function result = apr15_RAM_dependence_case(P,run_index,output_path)
%APR15_RAM_DEPENDENCE_CASE Run one fresh-process RAM measurement case.
%
% This helper is used by apr15_RAM_dependence.m so each particle count is
% measured in a clean external MATLAB process.

if nargin < 2
    error('apr15_RAM_dependence_case requires P and run_index.');
end

if nargin < 3
    output_path = '';
end

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));
set(0,'DefaultFigureVisible','off');

geom = struct();
geom.domain = 'boxed';
geom.phi = 0.55;
geom.rad = 1;
geom.min_gap = 1e-3;
geom.n_sweeps = 20;

solver = struct();
solver.N_c = 60;
solver.N_f = 150;
solver.N_peanut = 400;
solver.delta_pair = 0.2;
solver.gmres_tol = 1e-8;
solver.maxit = 1000;

geom_seed = 1000*P + run_index;

geom_opt = struct();
geom_opt.domain = geom.domain;
geom_opt.phi = geom.phi;
geom_opt.rad = geom.rad;
geom_opt.min_gap = geom.min_gap;
geom_opt.n_sweeps = geom.n_sweeps;
geom_opt.rng_seed = geom_seed;
geom_opt.visualise = false;

[q,geom_meta] = random_discs_mc(P,geom_opt);

opt = get2Dparams(P,solver.N_c,solver.N_f);
opt.rad = geom.rad;
opt.delta_pair = solver.delta_pair;
opt.N_peanut = solver.N_peanut;
opt.gmres_tol = solver.gmres_tol;
opt.maxit = solver.maxit;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.debug = 0;
opt.gmres_verbose = 0;
opt.surface_error_mode = 'rel';
opt.reuse_pair_basis_by_sep = false;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.get_bndry_field = 1;
opt.RAM_check = true;

F = zeros(P,2);
T = ones(P,1);

tic;
[~,sol] = solve_mob_peanut_enhanced(q,F,T,opt);
solve_time = toc;

ram = sol.ram_estimate;

result = struct();
result.P = P;
result.run_index = run_index;
result.geom_seed = geom_seed;
result.geom_meta = geom_meta;
result.it = sol.it;
result.gmres_unknowns = sol.gmres_unknowns;
result.rel_res = sol.rel_res;
result.solve_time = solve_time;
result.baseline_bytes = ram.baseline_bytes;
result.overall_peak_bytes = ram.overall_peak_bytes;
result.overall_delta_bytes = ram.overall_delta_bytes;
result.precomp_peak_bytes = ram.precomp.peak_bytes;
result.solve_peak_bytes = ram.solve.peak_bytes;
result.postprocess_peak_bytes = ram.postprocess.peak_bytes;
result.precomp_delta_bytes = ram.precomp.delta_bytes;
result.solve_delta_bytes = ram.solve.delta_bytes;
result.postprocess_delta_bytes = ram.postprocess.delta_bytes;

if ~isempty(output_path)
    save(output_path,'result');
end
end
