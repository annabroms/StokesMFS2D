function report = test_elastance_cmap_interpolation()
%TESTELASTANCECMAPINTERPOLATION Train/load and test one full-map model.
%
% Edit the configuration block for new parameter studies. If the
% parameter-keyed model is absent, interactive runs ask before training;
% batch runs never train implicitly.

%% Configuration
interpolation_tol = 1e-5;
volt_charge_interp_tol = 1e-7;
R = 2;
P = 3;
N_cmap = 80;
N_f = 150;
N_peanut = 400;
delta_min = 1e-3*R;
delta_pair = 0.2*R;
refit_cases = [false true];
refit_solver_N_c = 60;
gmres_tol = 1e-9;
solve_error_limit = 1e-5;
operator_audit_factor = 2;
rng_seed = 20260928;

opt = getLaplace2Dparams(P,R,N_cmap,N_f);
opt.N_cmap = N_cmap;
opt.N_peanut = N_peanut;
opt.delta_pair = delta_pair;
opt.smallest_delta = delta_min;
opt.use_interpolation = 'full';
opt.interpolation_tol = interpolation_tol;
opt.volt_charge_interp_tol = volt_charge_interp_tol;
opt.show_counter = false;
opt.parallel_precomp = false;
opt.get_bndry_field = true;
opt.use_fmm = false;
opt.gmres_tol = gmres_tol;
opt.gmres_verbose = 0;
opt.get_solve_time = false;
opt.visualise_sol = false;

fprintf('\n=== test_elastance_cmap_interpolation ===\n');
opt = prepareLaplaceCmapInterpolation(opt,'elastance');
model = opt.interpolation_model;
assert(strcmp(getLaplaceCmapModelProblem(model),'elastance'));
assert(strcmp(model.mode,'full'));

%% Independent canonical-map audit at endpoints and off-node panel points
audit_gaps = zeros(4*model.n_panels,1);
for k = 1:model.n_panels
    panel = model.panels(k);
    index = 4*(k-1)+(1:4);
    audit_gaps(index) = [panel.delta_lo; ...
        panel.delta_lo+0.271*(panel.delta_hi-panel.delta_lo); ...
        panel.delta_lo+0.733*(panel.delta_hi-panel.delta_lo); ...
        panel.delta_hi];
end
audit_gaps = unique(audit_gaps);
grids = struct('rbase_in_c',model.rbase_in_c, ...
    'rbase_in_f',model.rbase_in_f, ...
    'rout_base_f',model.rout_base_f);
exact_map_opt = removeModelFields(opt);
exact_map_opt.use_interpolation = 'none';
map_error = zeros(size(audit_gaps));
volt_charge_map_error = zeros(size(audit_gaps));
for k = 1:numel(audit_gaps)
    exact = buildCanonicalLaplacePairMap( ...
        audit_gaps(k),exact_map_opt,grids,'elastance');
    approximation = evaluateLaplaceCmapInterpolator( ...
        model,audit_gaps(k));
    map_error(k) = relativeSpectralError( ...
        approximation.Cmap,exact.Cmap);
    volt_charge_map_error(k) = relativeSpectralError( ...
        approximation.Cmap_QV,exact.Cmap_QV);
end
assert(max(map_error) <= operator_audit_factor*interpolation_tol, ...
    'Independent elastance Cmap audit exceeded its tolerance allowance.');
assert(max(volt_charge_map_error) <= ...
        operator_audit_factor*volt_charge_interp_tol, ...
    'Independent elastance Cmap_QV audit exceeded its tolerance allowance.');

%% Complete-solve comparisons
q = [0;(2*R+0.02*R)*exp(1i*0.37); ...
    (2*R+0.18*R)*exp(1i*2.4)];
rng(rng_seed,'twister');
Q_body = randn(P,1);
Q_body = Q_body-mean(Q_body);
case_template = struct('refit',[],'N_c',[],'voltage_error',[], ...
    'maxres_exact',[],'maxres_interpolated',[],'iterations',[]);
cases = repmat(case_template,numel(refit_cases),1);
for index = 1:numel(refit_cases)
    solve_opt = opt;
    solve_opt.refit = logical(refit_cases(index));
    if solve_opt.refit
        solve_opt.N_c = refit_solver_N_c;
    else
        solve_opt.N_c = N_cmap;
    end
    exact_opt = removeModelFields(solve_opt);
    exact_opt.use_interpolation = 'none';
    [v_exact,sol_exact] = solve_elast_peanut(q,Q_body,exact_opt);
    [v_interpolated,sol_interpolated] = ...
        solve_elast_peanut(q,Q_body,solve_opt);
    voltage_error = norm(v_interpolated-v_exact)/max(norm(v_exact),eps);
    assert(voltage_error <= solve_error_limit, ...
        'Interpolated elastance voltage error exceeded %.3e.', ...
        solve_error_limit);
    cases(index).refit = solve_opt.refit;
    cases(index).N_c = solve_opt.N_c;
    cases(index).voltage_error = voltage_error;
    cases(index).maxres_exact = sol_exact.maxres;
    cases(index).maxres_interpolated = sol_interpolated.maxres;
    cases(index).iterations = [sol_exact.it sol_interpolated.it];
end

%% Reproducible close-packed solve comparison
packed_P = 8;
packing_opt = struct('rad',R,'domain','boxed','phi',0.62, ...
    'min_gap',delta_min*(1+1e-6),'n_sweeps',120, ...
    'rng_seed',rng_seed+1,'visualise',false);
[q_packed,packed_geometry] = random_discs_mc(packed_P,packing_opt);
rng(rng_seed+1,'twister');
Q_packed = randn(packed_P,1);
Q_packed = Q_packed-mean(Q_packed);
packed_opt = opt;
packed_opt.P = packed_P;
packed_opt.N_c = N_cmap;
packed_opt.refit = false;
packed_opt.get_bndry_field = false;
packed_exact_opt = removeModelFields(packed_opt);
packed_exact_opt.use_interpolation = 'none';
[v_packed_exact,~] = solve_elast_peanut( ...
    q_packed,Q_packed,packed_exact_opt);
[v_packed_interpolated,~] = solve_elast_peanut( ...
    q_packed,Q_packed,packed_opt);
packed_voltage_error = norm(v_packed_interpolated-v_packed_exact)/ ...
    max(norm(v_packed_exact),eps);
assert(packed_voltage_error <= solve_error_limit, ...
    'Close-packed elastance voltage error exceeded %.3e.', ...
    solve_error_limit);

%% Precomputed big-sparse parity on the equal-grid case
standard_opt = opt;
standard_opt.N_c = N_cmap;
standard_opt.refit = false;
standard_opt.get_bndry_field = false;
[v_standard,~] = solve_elast_peanut(q,Q_body,standard_opt);
sparse_opt = standard_opt;
sparse_opt.use_big_sparse = true;
sparse_opt.lap_big_sparse_build_mode = 'precomputed';
[v_sparse,~] = solve_elast_peanut(q,Q_body,sparse_opt);
sparse_error = norm(v_sparse-v_standard)/max(norm(v_standard),eps);
assert(sparse_error <= 1e-10, ...
    'Precomputed big-sparse elastance interpolation disagrees with the standard path.');

%% Negative contract checks
bad_opt = opt;
bad_opt.use_interpolation = 'reduced';
assertThrows(@() configureLaplaceCmapInterpolation( ...
    bad_opt,'elastance'),'configureLaplaceCmapInterpolation:ElastanceFullOnly');
bad_model_opt = opt;
bad_model_opt.interpolation_model = model;
bad_model_opt.interpolation_model.problem = 'capacitance';
bad_model_opt.interpolation_model.kind = 'laplace_capacitance_cmap_alpha';
assertThrows(@() prepareLaplaceCmapInterpolation( ...
    bad_model_opt,'elastance'), ...
    'prepareLaplaceCmapInterpolation:IncompatibleModel');

changed_opt = opt;
changed_opt.interpolation_tol = opt.interpolation_tol/2;
assertThrows(@() prepareLaplaceCmapInterpolation( ...
    changed_opt,'elastance'), ...
    'prepareLaplaceCmapInterpolation:IncompatibleModel');
assertThrows(@() evaluateLaplaceCmapInterpolator( ...
    model,model.delta_max+max(1,model.delta_max)), ...
    'evaluateLaplaceCmapInterpolator:GapOutOfRange');

missing_model_opt = removeModelFields(opt);
assertThrows(@() solve_elast_peanut(q,Q_body,missing_model_opt), ...
    'solve_elast_peanut:ModelNotPrepared');
streaming_opt = opt;
streaming_opt.use_big_sparse = true;
streaming_opt.get_bndry_field = false;
streaming_opt.lap_big_sparse_build_mode = 'streaming';
assertThrows(@() solve_elast_peanut(q,Q_body,streaming_opt), ...
    'solve_elast_peanut:InterpolationSparseRequiresPrecomputed');

legacy_alias_opt = removeModelFields(opt);
legacy_alias_opt = rmfield(legacy_alias_opt,'volt_charge_interp_tol');
legacy_alias_opt.charge_interpolation_tol = volt_charge_interp_tol;
[legacy_alias_opt,~] = configureLaplaceCmapInterpolation( ...
    legacy_alias_opt,'elastance');
assert(legacy_alias_opt.volt_charge_interp_tol == ...
    volt_charge_interp_tol);
conflicting_alias_opt = opt;
conflicting_alias_opt.charge_interpolation_tol = ...
    2*opt.volt_charge_interp_tol;
assertThrows(@() configureLaplaceCmapInterpolation( ...
    conflicting_alias_opt,'elastance'), ...
    'configureLaplaceCmapInterpolation:ConflictingTolerance');

report = struct('configuration',struct( ...
    'problem','elastance','mode','full','P',P,'R',R, ...
    'interpolation_tol',interpolation_tol, ...
    'volt_charge_interp_tol',volt_charge_interp_tol, ...
    'N_cmap',N_cmap,'N_f',N_f,'N_peanut',N_peanut, ...
    'delta_min',delta_min,'delta_pair',delta_pair), ...
    'model',model,'q',q,'Q_body',Q_body,'audit_gaps',audit_gaps, ...
    'map_error',map_error, ...
    'volt_charge_map_error',volt_charge_map_error, ...
    'cases',cases,'packed_geometry',packed_geometry, ...
    'q_packed',q_packed,'Q_packed',Q_packed, ...
    'packed_voltage_error',packed_voltage_error, ...
    'sparse_error',sparse_error);
fprintf('Full elastance interpolation test completed successfully.\n');
end

function opt = removeModelFields(opt)
fields = {'interpolation_model','interpolation_model_file'};
present = fields(isfield(opt,fields));
if ~isempty(present)
    opt = rmfield(opt,present);
end
end

function value = relativeSpectralError(approximation,exact)
value = norm(approximation-exact,2)/max(norm(exact,2),eps);
end

function assertThrows(action,identifier)
thrown = false;
try
    action();
catch exception
    thrown = strcmp(exception.identifier,identifier);
end
assert(thrown,'Expected error %s.',identifier);
end
