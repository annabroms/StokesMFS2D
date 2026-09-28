function result = sep23_laplace_shifted_layer_packing(opts)
%SEP23_LAPLACE_SHIFTED_LAYER_PACKING Off-grid staggered-layer diagnostic.
%
% Test pair-basis reuse in a staggered, multi-layer packing whose layer
% shift makes the close-pair directions genuinely off the Fourier grids.
% This is not a rigidly rotated hexagonal lattice: the particle grids stay
% fixed, while each successive layer is translated by row_dx+i*row_dy.
%
% This complements apr26_reuse_sep_rotation_compare.m.  The Apr 26 script
% sweeps the layer shift with random boundary samples as an operator stress
% test.  Here the data are fixed smooth functions, their one-body and
% pair-output Fourier tails are measured, and the production capacitance
% solver is also compared with and without pair-basis reuse.
%
% The default inter-layer angle is chosen halfway between two N_c Fourier
% grid angles near 60 degrees.  Thus the horizontal pairs are grid-aligned,
% while the inter-layer pairs exercise fractional rotations in an otherwise
% realistic close packing.
%
% Syntax:
%   result = sep23_laplace_shifted_layer_packing()
%   result = sep23_laplace_shifted_layer_packing(opts)
%
% The optional opts struct can contain n_rows, n_cols, R, delta, N_c, N_f,
% N_peanut, offgrid_fraction, alpha, run_full_solve, plotfig,
% field_tolerance, and tail_tolerance.
%
% Anna Broms, Sep 23, 2026

if nargin < 1 || isempty(opts)
    opts = struct();
elseif ~isstruct(opts)
    error('opts must be a struct when supplied.');
end

script_name = mfilename;
script_date = 'Sep 23, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s (%s) ===\n',script_name,script_date);

% -------------------------------------------------------------------------
% Configuration.
% -------------------------------------------------------------------------
n_rows = getOptField(opts,'n_rows',3);
n_cols = getOptField(opts,'n_cols',4);
R = getOptField(opts,'R',2);
delta = getOptField(opts,'delta',1e-3);
N_c = getOptField(opts,'N_c',60);
N_f = getOptField(opts,'N_f',80);
N_peanut = getOptField(opts,'N_peanut',160);
offgrid_fraction = getOptField(opts,'offgrid_fraction',0.33);
alpha = getOptField(opts,'alpha',[]);
run_full_solve = logical(getOptField(opts,'run_full_solve',true));
plotfig = logical(getOptField(opts,'plotfig',true));
field_tolerance = getOptField(opts,'field_tolerance',1e-6);
tail_tolerance = getOptField(opts,'tail_tolerance',1e-8);

assert(mod(N_c,2)==0 && mod(N_f,2)==0, ...
    'Even N_c and N_f are required for the Nyquist diagnostic.');
assert(n_rows >= 2 && n_cols >= 2, ...
    'Use at least two rows and two columns.');

d = 2*R*(1+delta);
coarse_step = 2*pi/N_c;
if isempty(alpha)
    % Pick the angle (k+offgrid_fraction)*2*pi/N_c nearest 60 degrees.
    angle_index = round((pi/3)/coarse_step-offgrid_fraction);
    target_angle = (angle_index+offgrid_fraction)*coarse_step;
    alpha = cos(target_angle);
end
% assert(alpha > 0 && alpha <= 0.5, ...
%     ['The present diagnostic expects 0 < alpha <= 0.5 so that the ', ...
%      'same-column inter-layer separation equals d.']);

row_dx = alpha*d;
row_dy = d*sqrt(1-alpha^2);
q = make_shifted_layers(n_rows,n_cols,d,row_dx,row_dy);
q = q-mean(q);
P = numel(q);
inter_layer_angle = atan2(row_dy,row_dx);

% Smooth, nontrivial body voltages: the restriction of an affine harmonic
% function to the particle centres.  Nearby bodies have different values,
% so the solve still resolves strong close-gap fields.
packing_scale = max(abs(q-mean(q)));
direction = exp(-0.37i);
v_body = 0.2 + 0.8*real(direction*(q-mean(q)))/packing_scale;

opt_base = make_options(P,R,N_c,N_f,N_peanut);
[~,~,~,~,~,pairs] = getEnhancedGrid(q,opt_base);
pair_table = describe_pairs(q,pairs,N_c,N_f,R);

fprintf(['Packing: %d x %d layers/columns, P=%d, alpha=%.9f\n', ...
    '  layer translation = %.9g + %.9gi, angle = %.6f deg\n', ...
    '  angle/coarse step = %.6f (fraction %.6f), pairs=%d\n'], ...
    n_rows,n_cols,P,alpha,row_dx,row_dy,inter_layer_angle*180/pi, ...
    inter_layer_angle/coarse_step, ...
    fractional_part(inter_layer_angle/coarse_step),size(pairs,1));
fprintf('  off-grid close pairs: coarse %d/%d, fine %d/%d\n\n', ...
    nnz(pair_table.offgrid_coarse),height(pair_table), ...
    nnz(pair_table.offgrid_fine),height(pair_table));

config = struct();
config.n_rows = n_rows;
config.n_cols = n_cols;
config.P = P;
config.R = R;
config.relative_gap = delta;
config.N_c = N_c;
config.N_f = N_f;
config.N_peanut = N_peanut;
config.alpha = alpha;
config.row_dx = row_dx;
config.row_dy = row_dy;
config.inter_layer_angle_degrees = inter_layer_angle*180/pi;
config.offgrid_fraction = offgrid_fraction;
config.field_tolerance = field_tolerance;
config.tail_tolerance = tail_tolerance;
config.run_full_solve = run_full_solve;

% -------------------------------------------------------------------------
% Build direct, FFT-reused, and oversampled-FFT-reused pair maps.
% -------------------------------------------------------------------------
variant_defs = [ ...
    struct('name','direct','reuse',false,'rotation_mode','fft'), ...
    struct('name','fft','reuse',true,'rotation_mode','fft'), ...
    struct('name','oversampled_fft','reuse',true, ...
        'rotation_mode','oversampled_fft')];
contexts = cell(size(variant_defs));
setup_times = zeros(size(variant_defs));
for iv = 1:numel(variant_defs)
    opt = opt_base;
    opt.reuse_pair_basis_by_sep = variant_defs(iv).reuse;
    opt.rotation_mode = variant_defs(iv).rotation_mode;
    fprintf('Building %-15s pair representation ... ',variant_defs(iv).name);
    [contexts{iv},setup_times(iv)] = build_context(q,opt);
    fprintf('%.2f s (%d separation groups)\n',setup_times(iv), ...
        contexts{iv}.basis.pair_cache.n_groups);
end

% Independent, angularly offset boundary points.  These are not any of the
% collocation grids used to construct the one-body or pair maps.
ncheck = 3*N_c+1;
tcheck = (0:ncheck-1)'*(2*pi/ncheck) + 0.371*(2*pi/ncheck);
rcheck = reshape(q.' + R*exp(1i*tcheck),[],1);

data_names = {'low_mode','harmonic_trace'};
transform_cases = repmat(empty_transform_case(),numel(data_names),1);
transform_spectra = repmat(struct(),numel(data_names),1);
for idata = 1:numel(data_names)
    data_name = data_names{idata};
    tau = make_smooth_data(data_name,q,contexts{1}.rbase_out,R);
    outputs = cell(size(variant_defs));
    for iv = 1:numel(variant_defs)
        outputs{iv} = apply_transform(tau,contexts{iv},rcheck);
    end

    direct = outputs{1};
    fft_out = outputs{2};
    over_out = outputs{3};
    item = empty_transform_case();
    item.data_name = data_name;
    item.onebody_tail_max = direct.onebody_spectrum.tail_max;
    item.onebody_tail_median = direct.onebody_spectrum.tail_median;
    item.onebody_nyquist_max = direct.onebody_spectrum.nyquist_max;
    item.pair_output_tail_max = direct.pair_spectrum.tail_max;
    item.pair_output_tail_median = direct.pair_spectrum.tail_median;
    item.pair_output_nyquist_max = direct.pair_spectrum.nyquist_max;
    item.final_tail_max = direct.final_spectrum.tail_max;
    item.final_nyquist_max = direct.final_spectrum.nyquist_max;
    item.coeff_error_fft = relative_error(fft_out.lambda,direct.lambda);
    item.coeff_error_over = relative_error(over_out.lambda,direct.lambda);
    item.pair_coeff_error_fft = relative_error( ...
        fft_out.pair_lambda,direct.pair_lambda);
    item.pair_coeff_error_over = relative_error( ...
        over_out.pair_lambda,direct.pair_lambda);
    item.u_corr_error_fft = relative_error(fft_out.u_corr,direct.u_corr);
    item.u_corr_error_over = relative_error(over_out.u_corr,direct.u_corr);
    item.complete_field_error_fft = relative_error( ...
        fft_out.field,direct.field);
    item.complete_field_error_over = relative_error( ...
        over_out.field,direct.field);
    item.charge_map_error_fft = relative_error( ...
        fft_out.pair_qv,direct.pair_qv);
    item.charge_map_error_over = relative_error( ...
        over_out.pair_qv,direct.pair_qv);
    item.fft_over_coeff_difference = relative_error( ...
        fft_out.lambda,over_out.lambda);
    item.fft_over_field_difference = relative_error( ...
        fft_out.field,over_out.field);
    transform_cases(idata) = item;

    transform_spectra(idata).data_name = data_name;
    transform_spectra(idata).mode = direct.final_spectrum.mode;
    transform_spectra(idata).onebody_envelope = ...
        direct.onebody_spectrum.envelope;
    transform_spectra(idata).pair_output_envelope = ...
        direct.pair_spectrum.envelope;
    transform_spectra(idata).final_envelope = ...
        direct.final_spectrum.envelope;
end

fprintf('\nSmooth-data transform on independent boundary points\n');
fprintf('%-16s %11s %11s %12s %12s %12s %12s\n', ...
    'data','1B tail','pair tail','coeff fft','field fft', ...
    'field over','fft-over');
for k = 1:numel(transform_cases)
    c = transform_cases(k);
    fprintf('%-16s %11.3e %11.3e %12.3e %12.3e %12.3e %12.3e\n', ...
        c.data_name,c.onebody_tail_max,c.pair_output_tail_max, ...
        c.coeff_error_fft,c.complete_field_error_fft, ...
        c.complete_field_error_over,c.fft_over_field_difference);
end

% -------------------------------------------------------------------------
% Production capacitance solve with smooth body voltages.
% -------------------------------------------------------------------------
solve_result = struct('enabled',false);
if run_full_solve
    solve_result = run_solver_comparison(q,v_body,opt_base,variant_defs,N_c);
end

max_transform_field_error = max([transform_cases.complete_field_error_fft]);
max_fft_over_difference = max([transform_cases.fft_over_field_difference]);
max_onebody_tail = max([transform_cases.onebody_tail_max]);
max_pair_tail = max([transform_cases.pair_output_tail_max]);

acceptance = struct();
acceptance.transform_field_pass = ...
    max_transform_field_error <= field_tolerance;
acceptance.fft_matches_oversampled = ...
    max_fft_over_difference <= 1e3*eps;
acceptance.onebody_spectrally_resolved = max_onebody_tail <= tail_tolerance;
acceptance.pair_output_spectrally_resolved = max_pair_tail <= tail_tolerance;
if run_full_solve
    acceptance.solver_charge_pass = ...
        solve_result.comparison.q_error_fft <= field_tolerance;
    acceptance.solver_boundary_pass = ...
        solve_result.runs(2).maxres <= max([field_tolerance, ...
        10*opt_base.gmres_tol,10*solve_result.runs(1).maxres]);
else
    acceptance.solver_charge_pass = NaN;
    acceptance.solver_boundary_pass = NaN;
end

result = struct();
result.config = config;
result.geometry = struct('q',q,'pairs',pairs,'pair_table',pair_table, ...
    'v_body',v_body);
result.setup_times = array2table(setup_times(:).', ...
    'VariableNames',{variant_defs.name});
result.transform = struct('cases',struct2table(transform_cases), ...
    'spectra',transform_spectra);
result.solve = solve_result;
result.acceptance = acceptance;
result.interpretation = interpret_result(acceptance,max_onebody_tail, ...
    max_pair_tail,field_tolerance);

print_conclusion(result);
if plotfig
    plot_result(result);
end

end

% =========================================================================
% Local helpers
% =========================================================================

function q = make_shifted_layers(n_rows,n_cols,d,row_dx,row_dy)
q = zeros(n_rows*n_cols,1);
for ir = 0:n_rows-1
    for ic = 0:n_cols-1
        q(ir*n_cols+ic+1) = ic*d + ir*(row_dx+1i*row_dy);
    end
end
end

function opt = make_options(P,R,N_c,N_f,N_peanut)
opt = getLaplace2Dparams(P,R,N_c,N_f);
opt.delta_pair = 0.2*R;
opt.N_peanut = N_peanut;
opt.cmap = true;
opt.compress_cmap = false;
opt.project_charge = false;
opt.reuse_pair_basis_by_sep = true;
opt.rotation_mode = 'fft';
opt.rotation_oversample = 8;
opt.parallel_precomp = false;
opt.check_rotations = false;
opt.use_fmm = false;
opt.use_big_sparse = false;
opt.get_bndry_field = true;
opt.show_counter = false;
opt.visualise_sol = false;
opt.visualise_grid = false;
opt.gmres_verbose = 0;
opt.gmres_tol = 1e-10;
opt.single_threaded = false;
opt.get_precomp_time = true;
opt.get_solve_time = true;
end

function [ctx,elapsed] = build_context(q,opt)
N_c = opt.N_c;
N_f = opt.N_f;
R = opt.rad;
nout = ceil(opt.a_c*N_c);

t_c = (0:N_c-1)'*(2*pi/N_c);
t_f = (0:N_f-1)'*(2*pi/N_f);
t_out = (0:nout-1)'*(2*pi/nout);
t_out_f = (0:ceil(opt.a_f*N_f)-1)'* ...
    (2*pi/ceil(opt.a_f*N_f));
rbase_in_c = opt.Rp_c*exp(1i*t_c);
rbase_in_f = opt.Rp_f*exp(1i*t_f);
rbase_out = R*exp(1i*t_out);
rout_base_f = R*exp(1i*t_out_f);

rvec_in = reshape(q.'+rbase_in_c,[],1);
rvec_out = reshape(q.'+rbase_out,[],1);
[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

timer = tic;
[Upf,Ypf,DC,YC,Cmap,Cmap_QV,pair_cache] = getPairBasisLaplace( ...
    q,rbase_in_c,rbase_in_f,rout_base_f,rbase_out, ...
    rimage_vec,refine,pairs,opt);
[U,Y] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out,[0 nout]);
elapsed = toc(timer);

% Use field assignments so cell-valued pair data cannot make MATLAB's
% struct constructor expand this into a nonscalar struct array.
geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.opt = opt;
geom.rvec_out = rvec_out;
geom.rcheck = rvec_out;
geom.q = q;
geom.pairs = pairs;
geom.rimage_vec = rimage_vec;
geom.rvec_in = rvec_in;
geom.pair_cache = pair_cache;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Upf = Upf;
basis.Ypf = Ypf;
basis.DC_all = DC;
basis.YC_all = YC;
basis.Cmap = Cmap;
basis.Cmap_QV = Cmap_QV;
basis.pair_cache = pair_cache;
basis.Nii = lapSLPmat(rbase_in_c,rbase_out);

ctx = struct();
ctx.geom = geom;
ctx.basis = basis;
ctx.rbase_out = rbase_out;
end

function tau = make_smooth_data(name,q,rbase_out,R)
P = numel(q);
nout = numel(rbase_out);
t = angle(rbase_out/R);
z = q.'+rbase_out;
scale = max(abs(q-mean(q)))+R;
switch name
    case 'low_mode'
        centre_amplitude = 1 + 0.15*real(q.'/scale) - ...
            0.08*imag(q.'/scale);
        blocks = (cos(2*t)+0.30*sin(3*t)+0.08*cos(6*t))* ...
            centre_amplitude + 0.12*real(q.'/scale);
    case 'harmonic_trace'
        w = (z-mean(q))/scale;
        blocks = real(1 + 0.35*w + 0.08*w.^2 + 0.02*w.^3);
    otherwise
        error('Unknown smooth-data family "%s".',name);
end
assert(isequal(size(blocks),[nout P]));
tau = blocks(:);
end

function out = apply_transform(tau,ctx,rcheck)
geom = ctx.geom;
geom.rcheck = rcheck;
[lambda,lambda_self,~,~,u_corr,pair_qv] = ...
    transform_lap_peanut(tau,geom,ctx.basis);
pair_lambda = lambda-lambda_self;
field = lapSLPfield(geom.rvec_in,rcheck,lambda,false)+u_corr;
P = numel(geom.q);
N = geom.opt.N_c;

out = struct();
out.lambda = lambda;
out.lambda_self = lambda_self;
out.pair_lambda = pair_lambda;
out.u_corr = u_corr;
out.field = field;
out.pair_qv = pair_qv;
out.onebody_spectrum = spectrum_metrics(reshape(lambda_self,N,P));
out.pair_spectrum = spectrum_metrics(reshape(pair_lambda,N,P));
out.final_spectrum = spectrum_metrics(reshape(lambda,N,P));
end

function table_out = describe_pairs(q,pairs,N_c,N_f,R)
n = size(pairs,1);
pair_id = (1:n).';
i = pairs(:,1);
j = pairs(:,2);
z = q(j)-q(i);
separation = abs(z);
surface_gap = separation-2*R;
angle_degrees = mod(angle(z),2*pi)*180/pi;
coarse_coordinate = mod(angle(z),2*pi)/(2*pi/N_c);
fine_coordinate = mod(angle(z),2*pi)/(2*pi/N_f);
coarse_fraction = fractional_part(coarse_coordinate);
fine_fraction = fractional_part(fine_coordinate);
offgrid_coarse = distance_to_integer(coarse_coordinate) > 1e-10;
offgrid_fine = distance_to_integer(fine_coordinate) > 1e-10;
table_out = table(pair_id,i,j,separation,surface_gap,angle_degrees, ...
    coarse_fraction,fine_fraction,offgrid_coarse,offgrid_fine);
end

function s = spectrum_metrics(blocks)
N = size(blocks,1);
F = fft(blocks,[],1);
k = [0:N/2 -N/2+1:-1].';
energy = sum(abs(F).^2,1);
energy = max(energy,realmin);
tail = sum(abs(F(abs(k)>N/4,:)).^2,1)./energy;
nyquist = abs(F(N/2+1,:)).^2./energy;

envelope = zeros(N/2+1,1);
envelope(1) = max(abs(F(1,:)),[],'all');
for m = 1:N/2-1
    envelope(m+1) = max(abs(F([m+1 N-m+1],:)),[],'all');
end
envelope(end) = max(abs(F(N/2+1,:)),[],'all');
envelope = envelope/max(max(envelope),realmin);

s = struct('tail_max',max(tail),'tail_median',median(tail), ...
    'nyquist_max',max(nyquist),'nyquist_median',median(nyquist), ...
    'mode',(0:N/2).','envelope',envelope);
end

function solve_result = run_solver_comparison(q,v_body,opt_base,defs,N_c)
fprintf('\nProduction solve with smooth affine body voltages\n');
runs = repmat(struct('name','','Q',[],'lambda_proxy',[], ...
    'it',0,'maxres',NaN,'wall_time',NaN,'precomp_time',[], ...
    'solve_time',[],'pair_precomp_stats',[],'spectrum',[]),numel(defs),1);
for iv = 1:numel(defs)
    opt = opt_base;
    opt.reuse_pair_basis_by_sep = defs(iv).reuse;
    opt.rotation_mode = defs(iv).rotation_mode;
    timer = tic;
    [Q,sol] = solve_cap_peanut(q,v_body,opt);
    runs(iv).name = defs(iv).name;
    runs(iv).Q = Q;
    runs(iv).lambda_proxy = sol.lambda_proxy;
    runs(iv).it = sol.it;
    runs(iv).maxres = sol.maxres;
    runs(iv).wall_time = toc(timer);
    runs(iv).precomp_time = sol.precomp_time;
    runs(iv).solve_time = sol.solve_time;
    runs(iv).pair_precomp_stats = sol.pair_precomp_stats;
    runs(iv).spectrum = spectrum_metrics( ...
        reshape(sol.lambda_proxy,N_c,numel(q)));
    fprintf('  %-15s it=%3d, maxres=%.3e, wall=%.2f s, tail=%.3e\n', ...
        runs(iv).name,runs(iv).it,runs(iv).maxres, ...
        runs(iv).wall_time,runs(iv).spectrum.tail_max);
end

direct = runs(1);
fft_run = runs(2);
over_run = runs(3);
comparison = struct();
comparison.q_error_fft = relative_error(fft_run.Q,direct.Q);
comparison.q_error_over = relative_error(over_run.Q,direct.Q);
comparison.coeff_error_fft = relative_error( ...
    fft_run.lambda_proxy,direct.lambda_proxy);
comparison.coeff_error_over = relative_error( ...
    over_run.lambda_proxy,direct.lambda_proxy);
comparison.fft_over_q_difference = relative_error( ...
    fft_run.Q,over_run.Q);
comparison.fft_over_coeff_difference = relative_error( ...
    fft_run.lambda_proxy,over_run.lambda_proxy);
fprintf(['  relative Q error: fft %.3e, oversampled %.3e; ', ...
    'fft-over %.3e\n'],comparison.q_error_fft, ...
    comparison.q_error_over,comparison.fft_over_q_difference);
fprintf('  relative coefficient error: fft %.3e, oversampled %.3e\n', ...
    comparison.coeff_error_fft,comparison.coeff_error_over);

solve_result = struct('enabled',true,'runs',runs, ...
    'comparison',comparison);
end

function item = empty_transform_case()
item = struct('data_name','', ...
    'onebody_tail_max',NaN,'onebody_tail_median',NaN, ...
    'onebody_nyquist_max',NaN, ...
    'pair_output_tail_max',NaN,'pair_output_tail_median',NaN, ...
    'pair_output_nyquist_max',NaN, ...
    'final_tail_max',NaN,'final_nyquist_max',NaN, ...
    'coeff_error_fft',NaN,'coeff_error_over',NaN, ...
    'pair_coeff_error_fft',NaN,'pair_coeff_error_over',NaN, ...
    'u_corr_error_fft',NaN,'u_corr_error_over',NaN, ...
    'complete_field_error_fft',NaN, ...
    'complete_field_error_over',NaN, ...
    'charge_map_error_fft',NaN,'charge_map_error_over',NaN, ...
    'fft_over_coeff_difference',NaN, ...
    'fft_over_field_difference',NaN);
end

function message = interpret_result(a,onebody_tail,pair_tail,field_tol)
if a.transform_field_pass && ...
        (isnan(a.solver_charge_pass) || a.solver_charge_pass)
    if a.onebody_spectrally_resolved && a.pair_output_spectrally_resolved
        message = ['The smooth one-body and pair-output coefficients are ', ...
            'spectrally resolved, and reuse agrees with direct construction.'];
    else
        message = sprintf([ ...
            'Physical fields agree to the requested %.1e tolerance even ', ...
            'though the measured tails are one-body %.3e and pair %.3e. ', ...
            'Judge reuse by the field error, not coefficient appearance alone.'], ...
            field_tol,onebody_tail,pair_tail);
    end
elseif ~a.onebody_spectrally_resolved
    message = ['The one-body coefficients already have a substantial ', ...
        'unresolved tail; Fourier transfer is not a safe representation ', ...
        'at this N_c.'];
else
    message = ['The one-body input is spectrally resolved but the reused ', ...
        'physical field does not match the direct construction. Inspect ', ...
        'the pair-output transfer/discrete-map equivariance.'];
end
end

function print_conclusion(result)
a = result.acceptance;
fprintf('\nAcceptance summary\n');
fprintf('  transform field <= tolerance:       %s\n',pass_text(a.transform_field_pass));
fprintf('  FFT agrees with oversampled FFT:    %s\n', ...
    pass_text(a.fft_matches_oversampled));
fprintf('  one-body coefficients resolved:     %s\n', ...
    pass_text(a.onebody_spectrally_resolved));
fprintf('  pair-output coefficients resolved:  %s\n', ...
    pass_text(a.pair_output_spectrally_resolved));
if result.solve.enabled
    fprintf('  production-solve charges agree:     %s\n', ...
        pass_text(a.solver_charge_pass));
    fprintf('  reused boundary residual acceptable:%s\n', ...
        pass_text(a.solver_boundary_pass));
end
fprintf('  Interpretation: %s\n',result.interpretation);
end

function plot_result(result)
q = result.geometry.q;
pairs = result.geometry.pairs;
R = result.config.R;
T = result.geometry.pair_table;

figure('Name','Sep 23 shifted-layer packing','Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
hold on
t = linspace(0,2*pi,160);
for k = 1:numel(q)
    plot(real(q(k)+R*exp(1i*t)),imag(q(k)+R*exp(1i*t)), ...
        'k-','LineWidth',0.8);
end
for row = 1:size(pairs,1)
    z = q(pairs(row,:));
    if T.offgrid_coarse(row)
        colour = [0.85 0.2 0.1];
    else
        colour = [0.15 0.35 0.8];
    end
    plot(real(z),imag(z),'-','Color',colour,'LineWidth',1.2);
end
axis equal tight
grid on
xlabel('x'); ylabel('y');
title(sprintf('Shifted layers, \\alpha=%.6f',result.config.alpha));

nexttile;
hold on
spec = result.transform.spectra(end);
semilogy(spec.mode,spec.onebody_envelope,'LineWidth',1.5);
semilogy(spec.mode,spec.pair_output_envelope,'LineWidth',1.5);
semilogy(spec.mode,spec.final_envelope,'LineWidth',1.5);
yline(result.config.tail_tolerance,':k');
grid on
xlabel('|Fourier mode|'); ylabel('normalised maximum amplitude');
legend('one-body','pair output','final','tail threshold', ...
    'Location','southwest');
title('Harmonic-trace coefficient spectra');

nexttile;
cases = result.transform.cases;
values = [cases.complete_field_error_fft cases.complete_field_error_over];
bar(values);
set(gca,'YScale','log','XTick',1:height(cases), ...
    'XTickLabel',cases.data_name);
yline(result.config.field_tolerance,':k');
ylabel('relative independent-boundary field error');
legend('fft','oversampled fft','Location','best');
grid on
title('Smooth-data transform');

nexttile;
if result.solve.enabled
    Q = [result.solve.runs.Q];
    plot(1:size(Q,1),Q,'o-','LineWidth',1.1);
    xlabel('particle'); ylabel('charge');
    legend({result.solve.runs.name},'Location','best');
    grid on
    title('Production capacitance solve');
else
    axis off
    text(0.5,0.5,'Full solve disabled','HorizontalAlignment','center');
end
end

function value = fractional_part(x)
value = mod(x,1);
value(abs(value-1)<100*eps) = 0;
end

function d = distance_to_integer(x)
f = fractional_part(x);
d = min(f,1-f);
end

function e = relative_error(test,reference)
e = norm(test-reference)/max(norm(reference),realmin);
end

function value = pass_text(tf)
if isnan(tf)
    value = 'not run';
elseif tf
    value = 'PASS';
else
    value = 'FAIL';
end
end
