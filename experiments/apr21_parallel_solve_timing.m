% Compare serial and parallel solve matvec timing for peanut mobility.
%
% Defaults are intentionally modest. Override P, pool_size, or n_repeats in
% the workspace before running this script for heavier local experiments.

close all;
clear;

script_name = mfilename;
script_date = 'Apr 21, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root,'startup.m'));
end

fprintf('=== %s (%s) ===\n',script_name,script_date);

if ~exist('P','var') || isempty(P)
    P = 100;
end
if ~exist('pool_size','var') || isempty(pool_size)
    pool_size = 8;
end
if ~exist('n_repeats','var') || isempty(n_repeats)
    n_repeats = 1;
end
if ~exist('results_path','var') || isempty(results_path)
    results_path = fullfile(repo_root,'experiments', ...
        'apr21_parallel_solve_timing_results.mat');
end

if ~exist('solver','var') || ~isstruct(solver)
    solver = struct();
end
solver = set_default_field(solver,'N_c',60);
solver = set_default_field(solver,'N_f',150);
solver = set_default_field(solver,'N_peanut',400);
solver = set_default_field(solver,'delta_pair',0.2);
solver = set_default_field(solver,'gmres_tol',1e-8);
solver = set_default_field(solver,'maxit',1000);
solver = set_default_field(solver,'use_direct',true);
solver = set_default_field(solver,'parallel_solve_chunk_size',16);

rng(210421,'twister');
q = grow_cluster(P,1e-3,2,1,false,false,false);
[n_close,~,~] = count_close_pairs(q,solver.delta_pair,1);
fprintf('Geometry: P=%d, close pairs=%d, pool_size=%d\n',P,n_close,pool_size);

F = zeros(P,2);
T = ones(P,1);

ensure_pool(pool_size);

results = repmat(init_result(),n_repeats,2);
for irun = 1:n_repeats
    fprintf('\nRepeat %d/%d\n',irun,n_repeats);
    [UW_serial,sol_serial] = run_variant(q,F,T,solver,false);
    [UW_parallel,sol_parallel] = run_variant(q,F,T,solver,true);

    results(irun,1) = pack_result('serial',sol_serial);
    results(irun,2) = pack_result('parallel',sol_parallel);

    serial_s = get_solve_total(sol_serial);
    parallel_s = get_solve_total(sol_parallel);
    speedup = serial_s/max(parallel_s,eps);
    uw_rel = norm(UW_parallel-UW_serial,inf) / ...
        max(1,norm(UW_serial,inf));
    lambda_rel = norm(sol_parallel.lambda_c-sol_serial.lambda_c,inf) / ...
        max(1,norm(sol_serial.lambda_c,inf));
    resvec_rel = compare_resvec(sol_serial.resvec,sol_parallel.resvec);

    fprintf('  serial   solve %.3fs, it=%d\n',serial_s,sol_serial.it);
    fprintf('  parallel solve %.3fs, it=%d, speedup %.2fx\n', ...
        parallel_s,sol_parallel.it,speedup);
    fprintf('  rel diff: UW %.3e, lambda %.3e, resvec %.3e\n', ...
        uw_rel,lambda_rel,resvec_rel);
end

save(results_path,'script_name','script_date','P','pool_size','n_repeats', ...
    'solver','n_close','results');
fprintf('\nSaved results to:\n  %s\n',results_path);

function [UW,sol] = run_variant(q,F,T,solver,parallel_solve)
P = numel(q);
opt = get2Dparams(P,solver.N_c,solver.N_f);
opt.rad = 1;
opt.delta_pair = solver.delta_pair;
opt.N_peanut = solver.N_peanut;
opt.gmres_tol = solver.gmres_tol;
opt.maxit = solver.maxit;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.debug = 0;
opt.gmres_verbose = 0;
opt.surface_error_mode = 'rel';
opt.reuse_pair_basis_by_sep = true;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.get_bndry_field = 0;
opt.get_precomp_time = true;
opt.get_solve_time = true;
opt.parallel_precomp = 0;
opt.parallel_solve = parallel_solve;
opt.use_direct = solver.use_direct;
opt.parallel_solve_chunk_size = solver.parallel_solve_chunk_size;
opt.RAM_check = 0;

[UW,sol] = solve_mob_peanut_enhanced(q,F,T,opt);
end

function s = set_default_field(s,name,value)
if ~isfield(s,name) || isempty(s.(name))
    s.(name) = value;
end
end

function ensure_pool(pool_size)
if pool_size <= 0
    return
end
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers ~= pool_size
    delete(pool);
    pool = [];
end
if isempty(pool)
    parpool('local',pool_size);
end
end

function result = init_result()
result = struct('variant','','solve_time',nan,'it',nan,'rel_res',nan, ...
    'abs_res',nan,'final_est_res',nan,'gmres_unknowns',nan, ...
    'parallel_solve_stats',struct());
end

function result = pack_result(variant,sol)
result = init_result();
result.variant = variant;
result.solve_time = get_solve_total(sol);
result.it = sol.it;
result.rel_res = sol.rel_res;
result.abs_res = sol.abs_res;
if isfield(sol,'resvec') && ~isempty(sol.resvec)
    result.final_est_res = sol.resvec(end);
end
result.gmres_unknowns = sol.gmres_unknowns;
if isfield(sol,'parallel_solve_stats')
    result.parallel_solve_stats = sol.parallel_solve_stats;
end
end

function value = get_solve_total(sol)
value = NaN;
if isfield(sol,'solve_time') && isstruct(sol.solve_time) && ...
        isfield(sol.solve_time,'total')
    value = sol.solve_time.total;
end
end

function rel = compare_resvec(a,b)
n = min(numel(a),numel(b));
if n == 0
    rel = NaN;
    return
end
rel = norm(a(1:n)-b(1:n),inf)/max(1,norm(a(1:n),inf));
end
