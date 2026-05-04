close all;
clear;

script_name = mfilename;
script_date = 'Apr 14, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end
data_root = get_data_root(repo_root);

fprintf('=== Random Discs Peanut Mobility (Apr 14, 2026) ===\n'); 

N_threads = str2double(getenv('SLURM_CPUS_PER_TASK'));
N_workers = str2double(getenv('SLURM_NTASKS'));
if ~isnan(N_workers)
    fprintf('Number of cores %u \n',N_workers);
    fprintf('Number of threads per core %u \n', N_threads);
else
    N_workers = 4; % run locally
    N_threads = 2; 
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool(N_workers);  
    elseif pool.NumWorkers ~= N_workers
        delete(gcp('nocreate'));      
        pool = parpool(N_workers);
    end  
end

spmd; maxNumCompThreads(N_threads); end


% Geometry
geom.P = 100;
geom.rad = 1;
geom.domain = 'boxed';
geom.phi = 0.65;
geom.min_gap = 1e-3;

% Monte Carlo generation
mc.n_sweeps = 30;
mc.rng_seed = 1;
mc.show_generation_plot = false;

% Applied loads
loads.rng_seed = 11;

% Solver
solver.N_c = 80;
solver.N_f = 40; %was 60
solver.N_peanut = 200; %was 400
solver.Nclust = 80; %was 100
solver.delta_pair = 0.2;
solver.gmres_tol = 1e-7; %was 1e-8; 
solver.maxit = 1000;
solver.get_precomp_time = true;


% Solver opt overrides
opt_cfg = struct();
opt_cfg.get_bndry_field = 0;
opt_cfg.RAM_check = 0;
opt_cfg.mob_big_sparse_build_mode = 'precomputed';
opt_cfg.mob_big_sparse_build_mode = 'streaming';
opt_cfg.use_big_sparse = 1; %switch here!
opt_cfg.mob_sparse_map_coarse = 1;
opt_cfg.parallel_precomp = 1; 
opt_cfg.mob_big_sparse_chunk_pairs = '';
opt_cfg.parallel_precomp_chunk_pairs = '';
opt_cfg.parallel_big_sparse_build = 1; 
opt_cfg.solve_threads = '';



% Data I/O
% This experiment stores or loads raw solver data only.
io.load_results = load_saved_results_default();
io.save_results = ~io.load_results;

io.data_filename = build_results_filename(script_name, geom.P, geom.phi, solver.N_c);

if io.load_results && io.save_results
    error('apr14_random_discs_peanut_mobility:ConflictingIO', ...
        'Choose either io.load_results or io.save_results, not both.');
end

if io.save_results && ~isfolder(data_root)
    mkdir(data_root);
end

main_data_path = fullfile(data_root,io.data_filename);

if io.load_results
    if ~isfile(main_data_path)
        error('apr14_random_discs_peanut_mobility:MissingDataFile', ...
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

    opt = build_solver_opt(geom.P, solver, geom, opt_cfg);

    rng(loads.rng_seed,'twister');
    F = randn(geom.P,2);
    F = F - mean(F,1);
    T = randn(geom.P,1);
    % T = ones(geom.P,1);
    % F = zeros(geom.P,2);

    tic;
    [UW, sol] = solve_mob_peanut_enhanced(q, F, T, opt);
    solve_time = toc;

    [U, W] = unpack_UW(UW);
    body_forcing = sqrt(sum(F.^2,2) + (T/geom.rad).^2);
    body_speed = sqrt(sum(U.^2,2)) + geom.rad*abs(W);
    body_residual = sol.body_rel_res_max;
    gmres_residual = final_gmres_residual(sol);

    results = struct();
    results.script_name = script_name;
    results.script_date = script_date;
    results.geom = geom;
    results.mc = mc;
    results.loads = loads;
    results.solver = solver;
    results.geom_meta = geom_meta;
    results.q = q;
    results.F = F;
    results.T = T;
    results.UW = UW;
    results.sol = sol;
    results.solve_time = solve_time;
    results.body_forcing = body_forcing;
    results.body_speed = body_speed;
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
sol = results.sol;
solve_time = results.solve_time;
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
fprintf('  iterations=%d, unknowns=%d, gmres_res=%.3e, surf_rel=%.3e, time=%.2fs\n', ...
    sol.it, sol.gmres_unknowns, gmres_residual, sol.rel_res, solve_time);
if isfield(sol,'precomp_time') && isstruct(sol.precomp_time) && isfinite(sol.precomp_time.total)
    fprintf(['  precomp total=%.2fs, one-body=%.2fs, ', ...
        'two-body/peanut=%.2fs\n'], ...
        sol.precomp_time.total, ...
        sol.precomp_time.one_body, ...
        sol.precomp_time.two_body_or_peanut);
end

if opt.get_bndry_field
    max_body_residual = max(body_residual);
    fprintf('  max(body_rel_res_max)=%.3e\n', max_body_residual);
    fprintf('  |max(body_rel_res_max) - sol.rel_res| = %.3e\n', ...
        abs(max_body_residual - sol.rel_res));
end

fprintf('Solver timing: fmm/total %0.3f\n', sol.solve_time.fmm/sol.solve_time.total);

function value = final_gmres_residual(sol)
if isfield(sol,'resvec') && ~isempty(sol.resvec)
    value = sol.resvec(end);
else
    value = nan;
end
end

function [U, W] = unpack_UW(UW)
U = [UW(1:3:end), UW(2:3:end)];
W = UW(3:3:end);
end

function filename = build_results_filename(script_name, P, phi, N_c)
filename = sprintf('%s_P%d_phi%.3f_Nc%d_results.mat', ...
    script_name, P, phi, N_c);
end

function data_root = get_data_root(repo_root)
override = getenv('STOKESMFS2D_DATA_DIR');
if isempty(override)
    data_root = fullfile(repo_root,'data');
else
    data_root = override;
    fprintf('Using data directory override: %s\n', data_root);
end
end

function enabled = load_saved_results_default()
enabled = false;
end

function results = load_results_file(data_path)
loaded = load(data_path,'results');
if ~isfield(loaded,'results')
    error('apr14_random_discs_peanut_mobility:BadDataFile', ...
        'The data file %s does not contain a ''results'' struct.', data_path);
end
results = loaded.results;
end

function opt = build_solver_opt(P, solver, geom, overrides)
if nargin < 4
    overrides = struct();
end

opt = get2Dparams(P, solver.N_c, solver.N_f);
opt.rad = geom.rad;
opt.delta_pair = solver.delta_pair;
opt.Nclust = solver.Nclust; 
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
opt.get_precomp_time = solver.get_precomp_time;
opt.parallel_precomp = 1;
opt = apply_struct_overrides(opt, overrides);

end

function target = apply_struct_overrides(target, overrides)
if nargin < 2 || isempty(overrides)
    return
end

names = fieldnames(overrides);
for k = 1:numel(names)
    target.(names{k}) = overrides.(names{k});
end
end
