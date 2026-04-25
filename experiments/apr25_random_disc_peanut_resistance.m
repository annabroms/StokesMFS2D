clc;

script_name = mfilename;
script_date = 'Apr 25, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')))
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Random Discs Peanut Resistance (Apr 25, 2026) ===\n');

% Geometry
geom.P = 10;
geom.rad = 1;
geom.domain = 'boxed';
geom.phi = 0.65;
geom.min_gap = 1e-3;

% Monte Carlo generation
mc.n_sweeps = 50;
mc.rng_seed = 1;
mc.show_generation_plot = false;

% Applied loads (prescribed velocities for resistance problem)
loads.rng_seed = 11;

% Solver
solver.N_c = 120;
solver.N_f = 60;
% solver.N_c = 60;
% solver.N_f = 150;
solver.N_peanut = 400;
solver.delta_pair = 0.2;
solver.gmres_tol = 1e-8;
solver.maxit = 1000;
solver.get_precomp_time = true;

% Sparse compression mode
% Options: 'combined' or 'factored'
% 'combined': one sparse M_u_corr acting on [lambda_self; pair_tau]
% 'factored': build M_u_cross and M_u_peanut instead
sparse_u_corr_mode = 'combined';
res_smaller_ucorr = 1;  % 0: combined velocity correction map; 1: structured dense payload
fprintf('NOTE: Using sparse compression mode: %s\n', sparse_u_corr_mode);
fprintf('NOTE: Using res_smaller_ucorr: %d (0=combined velocity correction map, 1=structured dense payload)\n', res_smaller_ucorr);

% Data I/O
% Use io.save_results = true to run the solve and store raw data.
% Use io.load_results = true to load the saved data without rerunning.
run_test = 1; 
if run_test
    io.load_results = false;
    io.save_results = true;
else
    io.load_results = true;
    io.save_results = false;
end

io.data_filename = build_results_filename(script_name, geom.P, geom.phi, solver.N_c);

if io.load_results && io.save_results
    error('apr25_random_disc_peanut_resistance:ConflictingIO', ...
        'Choose either io.load_results or io.save_results, not both.');
end

main_data_path = fullfile(repo_root,'data',io.data_filename);

if io.load_results
    if ~isfile(main_data_path)
        error('apr25_random_disc_peanut_resistance:MissingDataFile', ...
            'Could not find saved results file: %s', main_data_path);
    end
    results = load_results_file(main_data_path);
    fprintf('Loaded results from %s\n', main_data_path);
else
    geom_opt = struct();
    geom_opt.domain = geom.domain;
    geom_opt.phi = geom.phi;
    geom_opt.rad = geom.rad;
    geom_opt.min_gap = geom.min_gap;
    geom_opt.n_sweeps = mc.n_sweeps;
    geom_opt.rng_seed = mc.rng_seed;
    geom_opt.visualise = mc.show_generation_plot;

    [q, geom_meta] = random_discs_mc(geom.P, geom_opt);

    opt = build_solver_opt(geom.P, solver, geom, sparse_u_corr_mode, res_smaller_ucorr);
    opt.get_bndry_field = 1;
    opt.RAM_check = 0; 
   % opt.Nclust = 40; 

    rng(loads.rng_seed,'twister');
    U = randn(geom.P,2);
    U = U - mean(U,1);
    W = randn(geom.P,1);
    W = ones(geom.P,1);
    U = zeros(geom.P,2);

    tic;
    [FT, sol] = solve_res_peanut_enhanced(q, U, W, opt);
    solve_time = toc;

    [F, T] = unpack_FT(FT);
    body_velocity = sqrt(sum(U.^2,2)) + geom.rad*abs(W);
    body_force = sqrt(sum(F.^2,2) + (T/geom.rad).^2);
    body_residual = sol.rel_res;
    gmres_residual = final_gmres_residual(sol);

    results = struct();
    results.script_name = script_name;
    results.script_date = script_date;
    results.geom = geom;
    results.mc = mc;
    results.loads = loads;
    results.solver = solver;
    results.sparse_u_corr_mode = sparse_u_corr_mode;
    results.res_smaller_ucorr = res_smaller_ucorr;
    results.geom_meta = geom_meta;
    results.q = q;
    results.U = U;
    results.W = W;
    results.FT = FT;
    results.sol = sol;
    results.solve_time = solve_time;
    results.body_velocity = body_velocity;
    results.body_force = body_force;
    results.body_residual = body_residual;
    results.gmres_residual = gmres_residual;
    results.opt = opt;

    if io.save_results
        save(main_data_path,'results');
        fprintf('Saved results to %s\n', main_data_path);
    end
end

geom_run = results.geom;
geom_meta = results.geom_meta;
q = results.q;
U = results.U;
W = results.W;
FT = results.FT;
sol = results.sol;
solve_time = results.solve_time;
body_velocity = results.body_velocity;
body_force = results.body_force;
body_residual = results.body_residual;
gmres_residual = results.gmres_residual;
opt = results.opt;

fprintf('Geometry:\n');
fprintf('  domain=%s, P=%d, phi_target=%.6f, phi=%.6f, L=%.6f\n', ...
    geom_meta.domain, geom_run.P, geom_meta.phi_target, geom_meta.phi, geom_meta.L);
fprintf('  min allowed gap=%.3e, measured min gap=%.3e\n', ...
    geom_meta.min_gap, geom_meta.min_surface_gap);

fprintf('Solver:\n');
fprintf('  N_c=%d, N_f=%d, N_peanut=%d, gmres_tol=%.1e\n', ...
    opt.N_c, opt.N_f, opt.N_peanut, opt.gmres_tol);
fprintf('  Sparse U-correction mode: %s\n', results.sparse_u_corr_mode);
fprintf('  res_smaller_ucorr: %d (0=combined velocity correction map, 1=structured dense payload)\n', opt.res_smaller_ucorr);
fprintf('  iterations=%d, unknowns=%d, gmres_res=%.3e, surf_rel=%.3e, time=%.2fs\n', ...
    sol.it, sol.gmres_unknowns, gmres_residual, sol.rel_res, solve_time);
if isfield(sol,'precomp_time') && isstruct(sol.precomp_time) && isfinite(sol.precomp_time.total)
    fprintf(['  precomp total=%.2fs, one-body=%.2fs, ', ...
        'two-body/peanut=%.2fs\n'], ...
        sol.precomp_time.total, ...
        sol.precomp_time.one_body, ...
        sol.precomp_time.two_body_or_peanut);
end

if isempty(body_residual)
    error('apr25_random_disc_peanut_resistance:MissingResiduals', ...
        'solve_res_peanut_enhanced did not return sol.body_rel_res_max.');
end

max_body_residual = max(body_residual);
fprintf('  max(body_rel_res_max)=%.3e\n', max_body_residual);
fprintf('  |max(body_rel_res_max) - sol.rel_res| = %.3e\n', ...
    abs(max_body_residual - sol.rel_res));

fprintf('Results saved to: %s\n', main_data_path);

function value = final_gmres_residual(sol)
if isfield(sol,'resvec') && ~isempty(sol.resvec)
    value = sol.resvec(end);
else
    value = nan;
end
end

function [F, T] = unpack_FT(FT)
F = [FT(1:3:end), FT(2:3:end)];
T = FT(3:3:end);
end

function filename = build_results_filename(script_name, P, phi, N_c)
filename = sprintf('%s_P%d_phi%.3f_Nc%d_results.mat', ...
    script_name, P, phi, N_c);
end

function results = load_results_file(data_path)
loaded = load(data_path,'results');
if ~isfield(loaded,'results')
    error('apr25_random_disc_peanut_resistance:BadDataFile', ...
        'The data file %s does not contain a ''results'' struct.', data_path);
end
results = loaded.results;
end

function save_results_file(data_path, results)
save(data_path,'results');
fprintf('Saved results to %s\n', data_path);
end

function opt = build_solver_opt(P, solver, geom, sparse_u_corr_mode, res_smaller_ucorr)
opt = get2Dparams(P, solver.N_c, solver.N_f);
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
opt.show_counter = 1;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.get_bndry_field = 0;
opt.get_precomp_time = solver.get_precomp_time;
opt.parallel_precomp = 1;
opt.res_big_sparse_u_corr_mode = sparse_u_corr_mode;
opt.res_smaller_ucorr = res_smaller_ucorr;
opt.use_big_sparse = 1;
opt.single_threaded = 1; 
end

