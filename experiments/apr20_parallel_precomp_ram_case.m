function result = apr20_parallel_precomp_ram_case(P,run_index,variant,pool_size,output_path)
%APR20_PARALLEL_PRECOMP_RAM_CASE Run one RAM benchmark case in fresh MATLAB.
%
% This helper is launched by apr20_parallel_precomp_ram_benchmark. It keeps
% the numerical problem fixed while varying the pair-precompute mode.

if nargin < 1 || isempty(P)
    P = 100;
end
if nargin < 2 || isempty(run_index)
    run_index = 1;
end
if nargin < 3 || isempty(variant)
    variant = 'parallel_streamed_auto_slim';
end
if nargin < 4 || isempty(pool_size)
    pool_size = 0;
end
if nargin < 5
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
opt.get_precomp_time = true;
opt.RAM_check = true;

variant = char(variant);
switch variant
    case 'serial_full'
        opt.parallel_precomp = false;
        opt.get_bndry_field = 1;
        pool_size = 0;
    case 'parallel_streamed_full'
        opt.parallel_precomp = true;
        opt.get_bndry_field = 1;
    case 'parallel_streamed_auto_slim'
        opt.parallel_precomp = true;
        opt.get_bndry_field = 0;
    case 'parallel_streamed_auto_slim_nodense'
        opt.parallel_precomp = true;
        opt.get_bndry_field = 0;
        opt.use_dense = 0;
    otherwise
        error('apr20_parallel_precomp_ram_case:UnknownVariant', ...
            'Unknown benchmark variant: %s',variant);
end

if opt.parallel_precomp
    ensure_parallel_pool(pool_size);
end

F = zeros(P,2);
T = ones(P,1);

solve_start_posix = current_posix_time();
tic;
[~,sol] = solve_mob_peanut_enhanced(q,F,T,opt);
solve_time = toc;
solve_end_posix = current_posix_time();

ram = sol.ram_estimate;
stats = sol.pair_precomp_stats;

result = struct();
result.P = P;
result.run_index = run_index;
result.variant = variant;
result.pool_size_requested = pool_size;
result.geom_seed = geom_seed;
result.geom_meta = geom_meta;
result.solve_time = solve_time;
result.solve_start_posix = solve_start_posix;
result.solve_end_posix = solve_end_posix;
result.precomp_time = sol.precomp_time;
result.pair_precomp_stats = stats;
result.client_baseline_bytes = ram.baseline_bytes;
result.client_overall_peak_bytes = ram.overall_peak_bytes;
result.client_precomp_peak_bytes = ram.precomp.peak_bytes;
result.client_solve_peak_bytes = ram.solve.peak_bytes;
result.client_postprocess_peak_bytes = ram.postprocess.peak_bytes;

if isfield(stats,'payload_mode')
    result.payload_mode = stats.payload_mode;
else
    result.payload_mode = '';
end
if isfield(stats,'parallel_backend')
    result.parallel_backend = stats.parallel_backend;
else
    result.parallel_backend = '';
end
if isfield(stats,'max_inflight')
    result.max_inflight = stats.max_inflight;
else
    result.max_inflight = NaN;
end

if ~isempty(output_path)
    save(output_path,'result');
end
end

function ensure_parallel_pool(pool_size)
if pool_size <= 0
    error('apr20_parallel_precomp_ram_case:InvalidPoolSize', ...
        'Parallel variants require a positive pool_size.');
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

function t = current_posix_time()
t = posixtime(datetime('now','TimeZone','UTC'));
end
