function report = test_capacitance_cmap_interpolation()
%TEST_CAPACITANCE_CMAP_INTERPOLATION Train/load and test one model setting.
%
% Edit only the configuration block below for routine experiments.  The
% selected interpolation model is loaded from its parameter-keyed MAT
% file.  If it is absent, an interactive run asks before training it.
% Batch runs never train implicitly.  For P>3 the complete-solve test uses
% a reproducible random close packing; P<=3 uses a small fixed geometry.
%
% Anna Broms, Sep 27, 2026

%% Configuration: model accuracy and canonical discretization
interpolation_mode = 'reduced_noconst'; % preferred interpolated mode. Choose 'reduced' or 'full'
interpolation_tol = 1e-6;       % final reconstructed C action error
charge_interpolation_tol = 1e-8; % independent full C_Q action error

R = 2;
P = 40;                         % P>3 selects a random close packing
N_cmap = 80;
N_f = 150;
N_peanut = 400;
delta_min = 1e-3*R;
delta_pair = 0.2*R;
interpolation_model_file = ''; % empty: parameter-keyed file under data/
panel_count_candidates = [1 2 4 8];
node_candidates = 3:2:17;
validation_nodes = 33;

%% Configuration: complete-solve comparison
refit_cases = [false true];
refit_solver_N_c = 60;         % used only when refit=true
gmres_tol = 1e-9;
solve_error_limit = 1e-5;      % diagnostic regression limit, not trainer target
operator_audit_factor = 2;     % allowance for independent off-node checks
packing_phi = 0.65;
packing_sweeps = 200;
rng_seed = 20260927;
show_packing = false;

%% Model options
validateattributes(P,{'numeric'},{'scalar','integer','>=',2},mfilename,'P');
opt = getLaplace2Dparams(P,R,N_cmap,N_f);
opt.N_cmap = N_cmap;
opt.N_peanut = N_peanut;
opt.delta_pair = delta_pair;
opt.smallest_delta = delta_min;
opt.use_interpolation = interpolation_mode;
opt.interpolation_tol = interpolation_tol;
opt.charge_interpolation_tol = charge_interpolation_tol;
opt.interpolation_model_file = interpolation_model_file;
opt.interpolation_panel_count_candidates = panel_count_candidates;
opt.interpolation_node_candidates = node_candidates;
opt.interpolation_validation_nodes = validation_nodes;
opt.show_counter = false;
opt.parallel_precomp = false;
opt.get_bndry_field = true;
opt.use_fmm = false;
opt.gmres_tol = gmres_tol;
opt.gmres_verbose = 0;
opt.get_solve_time = false;
opt.visualise_sol = false;

fprintf('\n=== test_capacitance_cmap_interpolation ===\n');
fprintf(['Requested mode=%s, P=%d, C tolerance %.3e, ', ...
    'C_Q tolerance %.3e.\n'],interpolation_mode,P, ...
    interpolation_tol,charge_interpolation_tol);
opt = prepareLaplaceCmapInterpolation(opt);
model = opt.interpolation_model;

%% Canonical off-node operator audit
audit_gaps = R*[0.011;0.137];
audit_gaps = min(max(audit_gaps,model.delta_min),model.delta_max);
grids = struct('rbase_in_c',model.rbase_in_c, ...
    'rbase_in_f',model.rbase_in_f, ...
    'rout_base_f',model.rout_base_f);
exact_map_opt = opt;
exact_map_opt.use_interpolation = 'none';
exact_map_opt = removeModelFields(exact_map_opt);
map_error = zeros(size(audit_gaps));
charge_map_error = zeros(size(audit_gaps));
for k = 1:numel(audit_gaps)
    exact = buildCanonicalLaplacePairMap( ...
        audit_gaps(k),exact_map_opt,grids);
    approximation = evaluateLaplaceCmapInterpolator( ...
        model,audit_gaps(k));
    panel = model.panels(approximation.panel_index);
    switch model.mode
        case 'reduced'
            C_approximation = panel.Cref+ ...
                panel.U*approximation.B*panel.V';
        case 'reduced_noconst'
            C_approximation = panel.U*approximation.B*panel.V';
        otherwise
            C_approximation = approximation.Cmap;
    end
    map_error(k) = relativeSpectralError( ...
        C_approximation,exact.Cmap);
    charge_map_error(k) = relativeSpectralError( ...
        approximation.Cmap_QV,exact.Cmap_QV);
end
fprintf('Off-node C errors:   %s\n',mat2str(map_error',4));
fprintf('Off-node C_Q errors: %s\n',mat2str(charge_map_error',4));
assert(max(map_error) <= operator_audit_factor*interpolation_tol, ...
    'Independent C-map audit exceeded %.1f times interpolation_tol.', ...
    operator_audit_factor);
assert(max(charge_map_error) <= ...
        operator_audit_factor*charge_interpolation_tol, ...
    ['Independent C_Q-map audit exceeded %.1f times ', ...
     'charge_interpolation_tol.'],operator_audit_factor);

%% Complete capacitance geometry
if P > 3
    packing_opt = struct('rad',R,'domain','boxed','phi',packing_phi, ...
        'min_gap',delta_min*(1+1e-6),'n_sweeps',packing_sweeps, ...
        'rng_seed',rng_seed,'visualise',show_packing);
    [q,geometry] = random_discs_mc(P,packing_opt);
    geometry_kind = 'random close packing';
else
    fixed = [0;(2*R+0.02*R)*exp(1i*0.37); ...
        (2*R+0.18*R)*exp(1i*2.4)];
    q = fixed(1:P);
    geometry = struct('min_surface_gap',minimumSurfaceGap(q,R), ...
        'phi',NaN,'rng_seed',rng_seed);
    geometry_kind = 'fixed off-axis cluster';
end
rng(rng_seed,'twister');
v_body = 2*rand(P,1)-1;
fprintf('\nGeometry: %s, P=%d, minimum surface gap %.6g.\n', ...
    geometry_kind,P,geometry.min_surface_gap);

%% Exact versus the one selected interpolation mode
case_template = struct('refit',[],'N_c',[],'Q_exact',[], ...
    'Q_interpolated',[],'charge_error',[],'maxres_exact',[], ...
    'maxres_interpolated',[],'iterations',[]);
cases = repmat(case_template,numel(refit_cases),1);
for index = 1:numel(refit_cases)
    refit = logical(refit_cases(index));
    solve_opt = opt;
    solve_opt.refit = refit;
    if refit
        solve_opt.N_c = refit_solver_N_c;
    else
        solve_opt.N_c = N_cmap;
    end

    exact_opt = solve_opt;
    exact_opt.use_interpolation = 'none';
    exact_opt = removeModelFields(exact_opt);
    [Q_exact,sol_exact] = solve_cap_peanut(q,v_body,exact_opt);
    [Q_interpolated,sol_interpolated] = ...
        solve_cap_peanut(q,v_body,solve_opt);
    charge_error = norm(Q_interpolated-Q_exact)/max(norm(Q_exact),eps);

    fprintf(['refit=%d, N_c=%d: relative Q error %.3e; ', ...
        'surface residual exact/interpolated %.3e / %.3e; ', ...
        'iterations %d / %d.\n'],refit,solve_opt.N_c,charge_error, ...
        sol_exact.maxres,sol_interpolated.maxres, ...
        sol_exact.it,sol_interpolated.it);
    assert(charge_error < solve_error_limit, ...
        ['Complete-solve Q error %.3e exceeds the configured diagnostic ', ...
         'limit %.3e. Tighten interpolation_tol or inspect the model ', ...
         'summary.'],charge_error,solve_error_limit);

    cases(index).refit = refit;
    cases(index).N_c = solve_opt.N_c;
    cases(index).Q_exact = Q_exact;
    cases(index).Q_interpolated = Q_interpolated;
    cases(index).charge_error = charge_error;
    cases(index).maxres_exact = sol_exact.maxres;
    cases(index).maxres_interpolated = sol_interpolated.maxres;
    cases(index).iterations = [sol_exact.it sol_interpolated.it];
end

report = struct('configuration',struct( ...
    'interpolation_mode',interpolation_mode,'P',P,'R',R, ...
    'interpolation_tol',interpolation_tol, ...
    'charge_interpolation_tol',charge_interpolation_tol, ...
    'N_cmap',N_cmap,'N_f',N_f,'N_peanut',N_peanut, ...
    'delta_min',delta_min,'delta_pair',delta_pair), ...
    'model',model,'geometry_kind',geometry_kind,'geometry',geometry, ...
    'q',q,'v_body',v_body,'audit_gaps',audit_gaps, ...
    'map_error',map_error,'charge_map_error',charge_map_error, ...
    'cases',cases);
fprintf('Selected %s model test completed successfully.\n',model.mode);
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

function gap = minimumSurfaceGap(q,R)
gap = Inf;
for i = 1:numel(q)-1
    gap = min(gap,min(abs(q(i)-q(i+1:end)))-2*R);
end
end
