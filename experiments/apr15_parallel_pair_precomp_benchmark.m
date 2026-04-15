close all;
clc;

% Checks speedup when precomputations are done with a parfor loop instead
% of serially.

script_name = mfilename;
script_date = 'Apr 15, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Parallel Pair Precomp Benchmark (%s) ===\n', script_date);

cfg = struct();
cfg.n_runs = 3;
cfg.warmup_runs = 1;
cfg.rel_tol = 1e-6;
cfg.gaps = 1e-3;

pool = gcp('nocreate');
if isempty(pool)
    fprintf('Opening local parallel pool...\n');
    pool = parpool('local');
else
    fprintf('Using existing pool with %d workers.\n', pool.NumWorkers);
end

specs = [build_stokes_spec(cfg.gaps), build_laplace_spec(cfg.gaps)];
reuse_modes = [false true];
results = cell(0,1);

for ss = 1:numel(specs)
        spec = specs(ss);
    for rr = 1:numel(reuse_modes)
        reuse_pair_basis = reuse_modes(rr);
        result = benchmark_mode(spec, reuse_pair_basis, cfg);
        results{end+1,1} = result; %#ok<SAGROW>
        print_result(result);
    end
end

fprintf('\n=== Summary Table ===\n');
fprintf(['%-24s %-9s %8s %8s %8s %12s %12s %12s %12s %12s\n'], ...
    'case','branch','pairs','groups','pool', ...
    'serial_pair','parallel_pair','speedup','serial_2b','parallel_2b');
for kk = 1:numel(results)
    result = results{kk};
    fprintf(['%-24s %-9s %8d %8d %8d %12.3f %12.3f %12.3f %12.3f %12.3f\n'], ...
        result.case_name, result.branch, result.n_pairs, result.n_groups, ...
        result.pool_size, result.serial_pair_basis_median, ...
        result.parallel_pair_basis_median, result.pair_basis_speedup, ...
        result.serial_two_body_median, result.parallel_two_body_median);
end

function spec = build_stokes_spec(gaps)
P = 32;
rad = 1;
q = grow_cluster(P,gaps,2);

opt = get2Dparams(P,60,150);
opt.rad = rad;
opt.delta_pair = 0.2;
opt.N_peanut = 400;
opt.gmres_tol = 1e-8;
opt.maxit = 600;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.debug = 0;
opt.gmres_verbose = 0;
opt.get_bndry_field = 1;
opt.get_precomp_time = true;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.parallel_precomp = false;

rng(11,'twister');
F = randn(P,2);
F = F - mean(F,1);
T = randn(P,1);

spec = struct();
spec.name = 'stokes_mob_peanut';
spec.solve_fun = @solve_mob_peanut_enhanced;
spec.opt = opt;
spec.inputs = {q, F, T};
spec.output_label = 'UW';
end

function spec = build_laplace_spec(gaps)
P = 24;
rad = 2;
q = grow_cluster(P,gaps,2,rad);

opt = getLaplace2Dparams(P,rad,40,80);
opt.delta_pair = max(gaps) + 0.05*rad;
opt.N_peanut = 200;
opt.gmres_tol = 1e-8;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.debug = 0;
opt.use_fmm = true;
opt.gmres_verbose = 0;
opt.get_bndry_field = 1;
opt.get_precomp_time = true;
opt.show_counter = 0;
opt.cmap = 1;
opt.parallel_precomp = false;

rng(17,'twister');
Q_body = randn(P,1);
Q_body = Q_body - mean(Q_body);

spec = struct();
spec.name = 'laplace_elast_peanut';
spec.solve_fun = @solve_elast_peanut;
spec.opt = opt;
spec.inputs = {q, Q_body};
spec.output_label = 'v_body';
end

function result = benchmark_mode(spec, reuse_pair_basis, cfg)
opt_serial = spec.opt;
opt_serial.reuse_pair_basis_by_sep = reuse_pair_basis;
opt_serial.parallel_precomp = false;

opt_parallel = opt_serial;
opt_parallel.parallel_precomp = true;

for kk = 1:cfg.warmup_runs
    run_solver(spec, opt_serial);
    run_solver(spec, opt_parallel);
end

serial_runs = cell(cfg.n_runs, 1);
parallel_runs = cell(cfg.n_runs, 1);
for kk = 1:cfg.n_runs
    serial_runs{kk} = run_solver(spec, opt_serial);
    parallel_runs{kk} = run_solver(spec, opt_parallel);
end

serial_last = serial_runs{end};
parallel_last = parallel_runs{end};
branch_expected = ternary(reuse_pair_basis, 'canonical_group', 'per_pair');
output_rel_diff = relerr_inf(serial_last.output, parallel_last.output);
serial_branch_match = strcmp(serial_last.sol.pair_precomp_stats.branch, branch_expected);
parallel_branch_match = strcmp(parallel_last.sol.pair_precomp_stats.branch, branch_expected);
parallel_used = parallel_last.sol.pair_precomp_stats.used_parallel;
output_within_tol = output_rel_diff <= cfg.rel_tol;

result = struct();
result.case_name = spec.name;
result.branch = branch_expected;
result.n_pairs = parallel_last.sol.pair_precomp_stats.n_pairs;
result.n_groups = parallel_last.sol.pair_precomp_stats.n_groups;
result.pool_size = parallel_last.sol.pair_precomp_stats.pool_size;
result.serial_pair_basis = cellfun(@(s) s.pair_basis, serial_runs);
result.parallel_pair_basis = cellfun(@(s) s.pair_basis, parallel_runs);
result.serial_two_body = cellfun(@(s) s.two_body, serial_runs);
result.parallel_two_body = cellfun(@(s) s.two_body, parallel_runs);
result.serial_pair_basis_median = median(result.serial_pair_basis);
result.parallel_pair_basis_median = median(result.parallel_pair_basis);
result.serial_two_body_median = median(result.serial_two_body);
result.parallel_two_body_median = median(result.parallel_two_body);
result.pair_basis_speedup = result.serial_pair_basis_median / ...
    max(result.parallel_pair_basis_median, eps);
result.two_body_speedup = result.serial_two_body_median / ...
    max(result.parallel_two_body_median, eps);
result.output_rel_diff = output_rel_diff;
result.serial_iters = serial_last.sol.it;
result.parallel_iters = parallel_last.sol.it;
result.serial_branch_actual = serial_last.sol.pair_precomp_stats.branch;
result.parallel_branch_actual = parallel_last.sol.pair_precomp_stats.branch;
result.serial_branch_match = serial_branch_match;
result.parallel_branch_match = parallel_branch_match;
result.parallel_used = parallel_used;
result.output_within_tol = output_within_tol;
result.rel_tol = cfg.rel_tol;
end

function run_result = run_solver(spec, opt)
args = spec.inputs;
tic;
[output, sol] = spec.solve_fun(args{:}, opt);
solve_time = toc;

run_result = struct();
run_result.output = output;
run_result.sol = sol;
run_result.solve_time = solve_time;
run_result.pair_basis = sol.precomp_time.pair_basis;
run_result.two_body = sol.precomp_time.two_body_or_peanut;
end


function print_result(result)
fprintf('\nCase: %s (%s)\n', result.case_name, result.branch);
fprintf('  n_pairs=%d, n_groups=%d, pool=%d\n', ...
    result.n_pairs, result.n_groups, result.pool_size);
fprintf('  pair_basis median: serial=%.3fs, parallel=%.3fs, speedup=%.3fx\n', ...
    result.serial_pair_basis_median, result.parallel_pair_basis_median, ...
    result.pair_basis_speedup);
fprintf('  two_body median  : serial=%.3fs, parallel=%.3fs, speedup=%.3fx\n', ...
    result.serial_two_body_median, result.parallel_two_body_median, ...
    result.two_body_speedup);
fprintf('  output rel diff  : %.3e\n', result.output_rel_diff);
fprintf('  iterations       : serial=%d, parallel=%d\n', ...
    result.serial_iters, result.parallel_iters);
fprintf(['  checks           : serial_branch=%s (%s), parallel_branch=%s (%s), ', ...
    'parallel_used=%s, output_tol=%s (tol=%.1e)\n'], ...
    pass_fail(result.serial_branch_match), result.serial_branch_actual, ...
    pass_fail(result.parallel_branch_match), result.parallel_branch_actual, ...
    pass_fail(result.parallel_used), pass_fail(result.output_within_tol), ...
    result.rel_tol);
end

function value = relerr_inf(a,b)
value = norm(a-b,inf) / max(1,norm(b,inf));
end

function value = ternary(cond, true_val, false_val)
if cond
    value = true_val;
else
    value = false_val;
end
end

function label = pass_fail(flag)
if flag
    label = 'OK';
else
    label = 'CHECK';
end
end
