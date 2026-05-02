function result = test_mobility_big_sparse_streaming
%TEST_MOBILITY_BIG_SPARSE_STREAMING Small tests for streamed Stokes big sparse.

fprintf('--- test_mobility_big_sparse_streaming ---\n');

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

q = (2.06*(0:3)).';
base_opt = make_small_options(numel(q));
result = struct();
result.ram_estimates = compare_ram_estimates(q,base_opt);

result.matrix_cases = [];
case_id = 0;
for direct_u_corr = [true false]
    for matrix_free = [true false]
        for sparse_map_coarse = [false true]
        case_id = case_id + 1;
        opt = base_opt;
        opt.big_sparse_direct_u_corr = direct_u_corr;
        opt.use_matrix_free_Lc_pair = matrix_free;
        opt.mob_sparse_map_coarse = sparse_map_coarse;
        result.matrix_cases(case_id).direct_u_corr = direct_u_corr;
        result.matrix_cases(case_id).use_matrix_free_Lc_pair = matrix_free;
        result.matrix_cases(case_id).sparse_map_coarse = sparse_map_coarse;
        result.matrix_cases(case_id).max_sparse_diff = ...
            compare_build_modes(q,opt);
        fprintf(['  matrix case %d: direct=%d matrix_free=%d ', ...
            'mob_sparse_map_coarse=%d max_sparse_diff=%.3e\n'], ...
            case_id,direct_u_corr,matrix_free,sparse_map_coarse, ...
            result.matrix_cases(case_id).max_sparse_diff);
        end
    end
end

result.solve_cases = [];
case_id = 0;
for direct_u_corr = [true false]
    for matrix_free = [true false]
        for sparse_map_coarse = [false true]
        case_id = case_id + 1;
        opt = base_opt;
        opt.big_sparse_direct_u_corr = direct_u_corr;
        opt.use_matrix_free_Lc_pair = matrix_free;
        opt.mob_sparse_map_coarse = sparse_map_coarse;
        [solve_rel_diff,it_pre,it_stream] = compare_solves(q,opt);
        result.solve_cases(case_id).direct_u_corr = direct_u_corr;
        result.solve_cases(case_id).use_matrix_free_Lc_pair = matrix_free;
        result.solve_cases(case_id).sparse_map_coarse = sparse_map_coarse;
        result.solve_cases(case_id).solve_rel_diff = solve_rel_diff;
        result.solve_cases(case_id).precomputed_iterations = it_pre;
        result.solve_cases(case_id).streaming_iterations = it_stream;
        fprintf(['  solve case %d: direct=%d matrix_free=%d ', ...
            'mob_sparse_map_coarse=%d rel_diff=%.3e\n'], ...
            case_id,direct_u_corr,matrix_free,sparse_map_coarse, ...
            solve_rel_diff);
        end
    end
end

result.parallel_build_cases = struct([]);
case_id = 0;
build_modes = {'precomputed','streaming'};
for mm = 1:numel(build_modes)
    build_mode = build_modes{mm};
    for direct_u_corr = [true false]
        for matrix_free = [true false]
            for sparse_map_coarse = [false true]
            case_id = case_id + 1;
            opt = base_opt;
            opt.big_sparse_direct_u_corr = direct_u_corr;
            opt.use_matrix_free_Lc_pair = matrix_free;
            opt.mob_sparse_map_coarse = sparse_map_coarse;
            case_result = compare_parallel_build(q,opt,build_mode);
            if case_id == 1
                result.parallel_build_cases = case_result;
            else
                result.parallel_build_cases(case_id) = case_result;
            end
            fprintf(['  parallel build case %d: mode=%s direct=%d ', ...
                'matrix_free=%d mob_sparse_map_coarse=%d ', ...
                'sparse=%.3e transform=%.3e matvec=%.3e solve=%.3e\n'], ...
                case_id,build_mode,direct_u_corr,matrix_free,...
                sparse_map_coarse,...
                result.parallel_build_cases(case_id).max_sparse_diff,...
                result.parallel_build_cases(case_id).transform_rel_diff,...
                result.parallel_build_cases(case_id).matvec_rel_diff,...
                result.parallel_build_cases(case_id).solve_rel_diff);
            end
        end
    end
end

assert(max([result.matrix_cases.max_sparse_diff]) < 5e-10, ...
    'Sparse matrices differ between precomputed and streaming builds.');
assert(max([result.solve_cases.solve_rel_diff]) < 5e-7, ...
    'Streaming and precomputed mobility solves differ.');
assert(max([result.parallel_build_cases.max_sparse_diff]) < 5e-10, ...
    'Parallel big sparse build changed the sparse matrices.');
assert(max([result.parallel_build_cases.transform_rel_diff]) < 5e-10, ...
    'Parallel big sparse build changed the mobility sparse transform.');
assert(max([result.parallel_build_cases.matvec_rel_diff]) < 5e-10, ...
    'Parallel big sparse build changed the mobility sparse matvec.');
assert(max([result.parallel_build_cases.solve_rel_diff]) < 5e-8, ...
    'Parallel big sparse build changed the mobility solve output.');
assert(max([result.parallel_build_cases.resvec_diff]) < 5e-10, ...
    'Parallel big sparse build changed the mobility GMRES history.');
assert(all([result.parallel_build_cases.it_serial] == ...
        [result.parallel_build_cases.it_parallel]), ...
    'Parallel big sparse build changed the mobility GMRES iteration count.');
assert(all([result.parallel_build_cases.used_parallel]), ...
    'Parallel big sparse build did not use the parallel backend.');
result.boundary_precomputed = compare_boundary_precomputed(q,base_opt);
result.boundary_parallel_precomp = compare_boundary_parallel_precomp(q,base_opt);
result.boundary_parallel_big_sparse = ...
    compare_boundary_parallel_big_sparse(q,base_opt);
result.chunk_normalization = check_chunk_normalization(q,base_opt);
assert(result.boundary_precomputed.UW_rel_diff < 5e-8, ...
    'Precomputed big sparse with boundary postprocessing changed UW.');
assert(isfinite(result.boundary_precomputed.big_sparse_rel_res), ...
    'Boundary postprocessing with precomputed big sparse failed.');
assert(result.boundary_parallel_precomp.UW_rel_diff < 1e-3, ...
    'Parallel precompute changed precomputed boundary UW too much.');
assert(result.boundary_parallel_precomp.it_serial == ...
        result.boundary_parallel_precomp.it_parallel, ...
    'Parallel precompute changed precomputed boundary iteration count.');
assert(result.boundary_parallel_big_sparse.UW_rel_diff < 5e-8, ...
    'Parallel big sparse build changed precomputed boundary UW.');
assert(result.boundary_parallel_big_sparse.rel_res_diff < 5e-10, ...
    'Parallel big sparse build changed precomputed boundary residuals.');
assert(result.boundary_parallel_big_sparse.resvec_diff < 5e-10, ...
    'Parallel big sparse build changed precomputed boundary GMRES history.');
assert(result.boundary_parallel_big_sparse.it_serial == ...
        result.boundary_parallel_big_sparse.it_parallel, ...
    'Parallel big sparse build changed precomputed boundary iteration count.');
assert(result.boundary_parallel_big_sparse.parallel_used, ...
    'Boundary precomputed solve did not use the parallel big sparse build.');
assert(result.ram_estimates.streaming_retained_pair_basis_bytes == 0, ...
    'Streaming should not retain pair-basis maps.');
assert(result.ram_estimates.precomputed_retained_pair_basis_bytes > 0, ...
    'Precomputed mode should report retained pair-basis maps.');
assert(result.ram_estimates.precomputed_solver_peak_bytes > ...
        result.ram_estimates.streaming_solver_peak_bytes, ...
    'Precomputed solver peak estimate should exceed streaming estimate.');

assert_streaming_errors(q,base_opt);

fprintf('  PASSED\n');
end

function opt = make_small_options(P)
opt = get2Dparams(P,12,18);
opt.a_c = 1.25;
opt.delta_pair = 0.12;
opt.N_peanut = 48;
opt.gmres_tol = 1e-7;
opt.maxit = 250;
opt.get_bndry_field = 0;
opt.get_precomp_time = true;
opt.get_solve_time = true;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.use_big_sparse = true;
opt.parallel_precomp_chunk_pairs = 2;
opt.mob_big_sparse_chunk_pairs = 2;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.use_fmm = false;
end

function max_diff = compare_build_modes(q,opt)
[geom,basis] = build_test_data(q,opt);

geom_pre = geom;
geom_pre.opt.mob_big_sparse_build_mode = 'precomputed';
[big_pre,~] = buildMobPeanutBigSparse(geom_pre,basis);

geom_stream = geom;
geom_stream.opt.mob_big_sparse_build_mode = 'streaming';
basis_stream = basis;
basis_stream.Upf = [];
basis_stream.Ypf = [];
basis_stream.DC_all = [];
basis_stream.YC_all = [];
basis_stream.Cmap = [];
basis_stream.Cmap_FU = [];
basis_stream.pair_cache = struct('enabled',false);
[big_stream,~] = buildMobPeanutBigSparse(geom_stream,basis_stream);

diffs = zeros(4,1);
names = {'M_source_corr','M_pair_nonp','M_pair_proj','M_rbm_corr', ...
    'M_u_corr','M_u_cross','M_u_peanut'};
for k = 1:numel(names)
    diffs(k) = compare_optional_sparse(big_stream,big_pre,names{k});
end
max_diff = max(diffs);
end

function case_result = compare_parallel_build(q,opt,build_mode)
ensure_pool(2);

opt_serial = opt;
opt_serial.mob_big_sparse_build_mode = build_mode;
opt_serial.parallel_precomp = false;
opt_serial.parallel_big_sparse_build = false;

opt_parallel = opt_serial;
opt_parallel.parallel_big_sparse_build = true;

[geom_serial,basis_serial] = build_test_data(q,opt_serial);
geom_parallel = geom_serial;
geom_parallel.opt = opt_parallel;
basis_parallel = basis_serial;

[basis_serial.big_sparse,stats_serial] = buildMobPeanutBigSparse( ...
    geom_serial,basis_serial);
[basis_parallel.big_sparse,stats_parallel] = buildMobPeanutBigSparse( ...
    geom_parallel,basis_parallel);

rng(13);
tau = randn(2*numel(geom_serial.rvec_out),1);
[lam_serial_x,lam_serial_y,lambda_self_serial,u_corr_serial] = ...
    transform_mob_peanut_big_sparse_stokes(tau,geom_serial,basis_serial);
[lam_parallel_x,lam_parallel_y,lambda_self_parallel,u_corr_parallel] = ...
    transform_mob_peanut_big_sparse_stokes(tau,geom_parallel,basis_parallel);

y_serial = matvec_peanut_big_sparse(tau,geom_serial,basis_serial);
y_parallel = matvec_peanut_big_sparse(tau,geom_parallel,basis_parallel);

[F,T] = build_test_loads();
[UW_serial,sol_serial] = solve_mob_peanut_enhanced(q,F,T,opt_serial);
[UW_parallel,sol_parallel] = solve_mob_peanut_enhanced(q,F,T,opt_parallel);

case_result = struct();
case_result.build_mode = build_mode;
case_result.direct_u_corr = opt.big_sparse_direct_u_corr;
case_result.use_matrix_free_Lc_pair = opt.use_matrix_free_Lc_pair;
case_result.sparse_map_coarse = opt.mob_sparse_map_coarse;
case_result.max_sparse_diff = compare_big_sparse_struct( ...
    basis_parallel.big_sparse,basis_serial.big_sparse);
case_result.transform_rel_diff = max([ ...
    relerr(lam_parallel_x,lam_serial_x), ...
    relerr(lam_parallel_y,lam_serial_y), ...
    relerr(lambda_self_parallel,lambda_self_serial), ...
    relerr(u_corr_parallel,u_corr_serial)]);
case_result.matvec_rel_diff = relerr(y_parallel,y_serial);
case_result.solve_rel_diff = relerr(UW_parallel,UW_serial);
case_result.resvec_diff = compare_vec_inf(sol_parallel.resvec,sol_serial.resvec);
case_result.it_serial = sol_serial.it;
case_result.it_parallel = sol_parallel.it;
case_result.used_parallel = stats_parallel.used_parallel;
case_result.parallel_backend = stats_parallel.parallel_backend;
case_result.chunk_pairs = stats_parallel.chunk_pairs;
case_result.n_tasks = stats_parallel.n_tasks;
end

function ram = compare_ram_estimates(q,opt)
n_pairs = count_close_pairs(q,opt.delta_pair,1);
N_check = ceil(opt.a_c*opt.N_c);

opt_pre = opt;
opt_pre.mob_big_sparse_build_mode = 'precomputed';
est_pre = estimateMobPeanutBigSparseRam( ...
    numel(q),opt.N_c,N_check,n_pairs,opt_pre);

opt_pre_boundary = opt_pre;
opt_pre_boundary.get_bndry_field = 1;
est_pre_boundary = estimateMobPeanutBigSparseRam( ...
    numel(q),opt.N_c,N_check,n_pairs,opt_pre_boundary);
assert(est_pre_boundary.retained_pair_basis_bytes > ...
        est_pre.retained_pair_basis_bytes, ...
    'Boundary precomputed estimate should include retained UB/YB pair factors.');

opt_stream = opt;
opt_stream.mob_big_sparse_build_mode = 'streaming';
est_stream = estimateMobPeanutBigSparseRam( ...
    numel(q),opt.N_c,N_check,n_pairs,opt_stream);

ram = struct();
ram.precomputed_retained_pair_basis_bytes = ...
    est_pre.retained_pair_basis_bytes;
ram.streaming_retained_pair_basis_bytes = ...
    est_stream.retained_pair_basis_bytes;
ram.precomputed_solver_peak_bytes = ...
    est_pre.solver_precompute_peak_bytes;
ram.streaming_solver_peak_bytes = ...
    est_stream.solver_precompute_peak_bytes;
ram.big_sparse_peak_diff_bytes = ...
    est_pre.big_sparse_peak_bytes - est_stream.big_sparse_peak_bytes;

fprintf(['  RAM estimate: retained pair basis precomputed %.2f MB, ', ...
    'streaming %.2f MB, solver peak diff %.2f MB\n'], ...
    bytes_to_mb(est_pre.retained_pair_basis_bytes), ...
    bytes_to_mb(est_stream.retained_pair_basis_bytes), ...
    (ram.precomputed_solver_peak_bytes - ...
     ram.streaming_solver_peak_bytes)/1024^2);
end

function mb = bytes_to_mb(bytes)
mb = bytes/1024^2;
end

function [solve_rel_diff,it_pre,it_stream] = compare_solves(q,opt)
[F,T] = build_test_loads();

opt_pre = opt;
opt_pre.mob_big_sparse_build_mode = 'precomputed';
[UW_pre,sol_pre] = solve_mob_peanut_enhanced(q,F,T,opt_pre);

opt_stream = opt;
opt_stream.mob_big_sparse_build_mode = 'streaming';
[UW_stream,sol_stream] = solve_mob_peanut_enhanced(q,F,T,opt_stream);

solve_rel_diff = norm(UW_stream-UW_pre,inf)/max(1,norm(UW_pre,inf));
it_pre = sol_pre.it;
it_stream = sol_stream.it;
end

function diff_value = compare_optional_sparse(a,b,field_name)
has_a = isfield(a,field_name);
has_b = isfield(b,field_name);
assert(has_a == has_b, ...
    'Streaming/precomputed sparse builds disagree on field %s.',field_name);
if ~has_a
    diff_value = 0;
    return
end
diff_value = norm(full(a.(field_name)-b.(field_name)),inf);
end

function result = compare_boundary_precomputed(q,opt)
[F,T] = build_test_loads();

opt_solve_grid = opt;
opt_solve_grid.use_big_sparse = true;
opt_solve_grid.mob_big_sparse_build_mode = 'precomputed';
opt_solve_grid.get_bndry_field = 0;
[UW_solve_grid,~] = solve_mob_peanut_enhanced(q,F,T,opt_solve_grid);

opt_boundary = opt;
opt_boundary.use_big_sparse = true;
opt_boundary.mob_big_sparse_build_mode = 'precomputed';
opt_boundary.get_bndry_field = 1;
[UW_boundary,sol_boundary] = solve_mob_peanut_enhanced(q,F,T,opt_boundary);

result = struct();
result.UW_rel_diff = norm(UW_boundary-UW_solve_grid,inf)/ ...
    max(1,norm(UW_solve_grid,inf));
result.big_sparse_rel_res = sol_boundary.rel_res;
fprintf(['  boundary precomputed: UW rel diff %.3e, ', ...
    'boundary rel_res %.3e\n'], ...
    result.UW_rel_diff,result.big_sparse_rel_res);
end

function result = compare_boundary_parallel_precomp(q,opt)
ensure_pool(2);

[F,T] = build_test_loads();

opt_serial = opt;
opt_serial.use_big_sparse = true;
opt_serial.mob_big_sparse_build_mode = 'precomputed';
opt_serial.get_bndry_field = 1;
opt_serial.parallel_precomp = false;
[UW_serial,sol_serial] = solve_mob_peanut_enhanced(q,F,T,opt_serial);

opt_parallel = opt_serial;
opt_parallel.parallel_precomp = true;
opt_parallel.parallel_precomp_chunk_pairs = 2;
[UW_parallel,sol_parallel] = solve_mob_peanut_enhanced(q,F,T,opt_parallel);

result = struct();
result.UW_rel_diff = norm(UW_parallel-UW_serial,inf)/max(1,norm(UW_serial,inf));
result.rel_res_diff = abs(sol_parallel.rel_res-sol_serial.rel_res);
result.resvec_diff = compare_vec_inf(sol_parallel.resvec,sol_serial.resvec);
result.it_serial = sol_serial.it;
result.it_parallel = sol_parallel.it;

fprintf(['  boundary parallel precomp: UW rel diff %.3e, ', ...
    'rel_res diff %.3e, resvec diff %.3e, it %d/%d\n'], ...
    result.UW_rel_diff,result.rel_res_diff,result.resvec_diff, ...
    result.it_serial,result.it_parallel);
end

function result = compare_boundary_parallel_big_sparse(q,opt)
ensure_pool(2);

[F,T] = build_test_loads();

opt_serial = opt;
opt_serial.use_big_sparse = true;
opt_serial.mob_big_sparse_build_mode = 'precomputed';
opt_serial.get_bndry_field = 1;
opt_serial.parallel_precomp = false;
opt_serial.parallel_big_sparse_build = false;
[UW_serial,sol_serial] = solve_mob_peanut_enhanced(q,F,T,opt_serial);

opt_parallel = opt_serial;
opt_parallel.parallel_big_sparse_build = true;
[UW_parallel,sol_parallel] = solve_mob_peanut_enhanced(q,F,T,opt_parallel);

result = struct();
result.UW_rel_diff = norm(UW_parallel-UW_serial,inf)/max(1,norm(UW_serial,inf));
result.rel_res_diff = abs(sol_parallel.rel_res-sol_serial.rel_res);
result.resvec_diff = compare_vec_inf(sol_parallel.resvec,sol_serial.resvec);
result.it_serial = sol_serial.it;
result.it_parallel = sol_parallel.it;
result.parallel_used = sol_parallel.big_sparse_stats.used_parallel;

fprintf(['  boundary parallel big sparse: UW rel diff %.3e, ', ...
    'rel_res diff %.3e, resvec diff %.3e, it %d/%d, used_parallel=%d\n'], ...
    result.UW_rel_diff,result.rel_res_diff,result.resvec_diff, ...
    result.it_serial,result.it_parallel,result.parallel_used);
end

function result = check_chunk_normalization(q,opt)
[~,~,~,~,~,pairs] = getEnhancedGrid(q,opt);
n_pairs = size(pairs,1);

parallel_only = opt;
parallel_only.parallel_precomp_chunk_pairs = 3;
parallel_only.mob_big_sparse_chunk_pairs = [];
[parallel_chunk,parallel_norm] = resolveMobilityChunkPairs( ...
    n_pairs,parallel_only,2);

mob_only = opt;
mob_only.parallel_precomp_chunk_pairs = [];
mob_only.mob_big_sparse_chunk_pairs = 5;
[mob_chunk,mob_norm] = resolveMobilityChunkPairs(n_pairs,mob_only,2);

both_equal = opt;
both_equal.parallel_precomp_chunk_pairs = 4;
both_equal.mob_big_sparse_chunk_pairs = 4;
[equal_chunk,equal_norm] = resolveMobilityChunkPairs(n_pairs,both_equal,2);

both_different = opt;
both_different.parallel_precomp_chunk_pairs = 2;
both_different.mob_big_sparse_chunk_pairs = 3;
did_error = false;
try
    resolveMobilityChunkPairs(n_pairs,both_different,2);
catch ME
    did_error = strcmp(ME.identifier, ...
        'resolveMobilityChunkPairs:ConflictingChunkPairs');
end

assert(parallel_chunk == 3 && parallel_norm.mob_big_sparse_chunk_pairs == 3, ...
    'Parallel-only chunk normalization failed.');
assert(mob_chunk == min(5,n_pairs) && ...
    mob_norm.parallel_precomp_chunk_pairs == min(5,n_pairs), ...
    'Big-sparse-only chunk normalization failed.');
assert(equal_chunk == min(4,n_pairs) && ...
    equal_norm.parallel_precomp_chunk_pairs == min(4,n_pairs) && ...
    equal_norm.mob_big_sparse_chunk_pairs == min(4,n_pairs), ...
    'Equal chunk normalization failed.');
assert(did_error,'Conflicting chunk settings should error.');

result = struct();
result.parallel_only_chunk = parallel_chunk;
result.mob_only_chunk = mob_chunk;
result.equal_chunk = equal_chunk;
result.conflict_error = did_error;

fprintf(['  chunk normalization: parallel-only=%d mob-only=%d ', ...
    'both=%d conflict=%d\n'], ...
    result.parallel_only_chunk,result.mob_only_chunk,...
    result.equal_chunk,result.conflict_error);
end

function [F,T] = build_test_loads()
F = [1 -0.5; -0.25 0.7; 0.4 -0.2; -1.15 0.0];
F = F - mean(F,1);
T = [0.2; -0.1; 0.05; -0.15];
end

function max_diff = compare_big_sparse_struct(a,b)
max_diff = 0;
names = {'M_source_corr','M_pair_nonp','M_pair_proj','M_rbm_corr', ...
    'M_u_corr','M_u_cross','M_u_peanut'};
for k = 1:numel(names)
    max_diff = max(max_diff,compare_optional_sparse(a,b,names{k}));
end
end

function [geom,basis] = build_test_data(q,opt)
q = q(:);
P = numel(q);
rad = 1;
N_c = opt.N_c;
N_f = opt.N_f;
a_c = opt.a_c;
Rp_c = opt.Rp_c;
Rp_f = opt.Rp_f;

tout_c_all = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c_all(1:end-1)';
rbase_out_c = rad*(cos(tout_c)+1i*sin(tout_c));
tin_c = linspace(0,2*pi,N_c+1)';
tin_c = tin_c(1:end-1);
rbase_in_c = Rp_c*cos(tin_c) + 1i*Rp_c*sin(tin_c);
tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
rbase_in_f = Rp_f*cos(tin_f) + 1i*Rp_f*sin(tin_f);

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

rvec_in_c = zeros(P*N_c,1);
rout = zeros(P*numel(rbase_out_c),1);
for k = 1:P
    rvec_in_c((k-1)*N_c+1:k*N_c) = q(k) + rbase_in_c;
    out_idx = (k-1)*numel(rbase_out_c)+1:k*numel(rbase_out_c);
    rout(out_idx) = q(k) + rbase_out_c;
end

svd_opts = struct('column_weight',false,'left_weight',false);
[U,Y,Lc] = getSelfPseudoMobilityStokes( ...
    1,q,rbase_in_c,rbase_out_c,[],[0,ceil(a_c*N_c)],svd_opts);

opt_pair = opt;
opt_pair.project_force = true;
opt_pair.project = true;
opt_pair.pair_basis_debug = 0;
opt_pair.rad = ones(P,1);
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU,pair_cache] = ...
    getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine, ...
    pairs,opt_pair,Lc{1},rbase_out_c,svd_opts);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.rvec_in = rvec_in_c;
geom.rimage_vec = rimage_vec;
geom.refine = refine;
geom.opt = opt_pair;
geom.opt.resistance = 0;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Lc = Lc{1};
basis.Upf = UB_all;
basis.Ypf = YB_all;
basis.DC_all = UC_all;
basis.YC_all = YC_all;
basis.Cmap = Cmap;
basis.Cmap_FU = Cmap_FU;
basis.Lc_pair = getILpair(Lc{1});
basis.pair_cache = pair_cache;
basis.Nii = stokSLPmat(rbase_in_c,rbase_out_c,1);
end

function assert_streaming_errors(q,opt)
F = [1 -0.5; -0.25 0.7; 0.4 -0.2; -1.15 0.0];
F = F - mean(F,1);
T = [0.2; -0.1; 0.05; -0.15];

opt_boundary = opt;
opt_boundary.mob_big_sparse_build_mode = 'streaming';
opt_boundary.get_bndry_field = 1;
did_error = false;
try
    solve_mob_peanut_enhanced(q,F,T,opt_boundary);
catch ME
    did_error = strcmp(ME.identifier, ...
        'solve_mob_peanut_enhanced:StreamingBoundaryUnsupported');
end
assert(did_error,'Streaming with get_bndry_field=1 should error.');

opt_reuse = opt;
opt_reuse.mob_big_sparse_build_mode = 'streaming';
opt_reuse.reuse_pair_basis_by_sep = true;
did_error = false;
try
    solve_mob_peanut_enhanced(q,F,T,opt_reuse);
catch ME
    did_error = strcmp(ME.identifier, ...
        'solve_mob_peanut_enhanced:StreamingPairReuseUnsupported');
end
assert(did_error,'Streaming with reuse_pair_basis_by_sep=1 should error.');
end

function value = relerr(a,b)
value = norm(a(:)-b(:),inf) / max(1,norm(b(:),inf));
end

function value = compare_vec_inf(a,b)
a = a(:);
b = b(:);
assert(numel(a) == numel(b), ...
    'Residual history length mismatch in mobility streaming comparison.');
value = norm(a-b,inf);
end

function ensure_pool(pool_size)
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers ~= pool_size
    delete(pool);
    pool = [];
end
if isempty(pool)
    parpool('local',pool_size);
end
end
