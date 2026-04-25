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

assert(max([result.matrix_cases.max_sparse_diff]) < 5e-10, ...
    'Sparse matrices differ between precomputed and streaming builds.');
assert(max([result.solve_cases.solve_rel_diff]) < 5e-8, ...
    'Streaming and precomputed mobility solves differ.');
result.boundary_precomputed = compare_boundary_precomputed(q,base_opt);
assert(result.boundary_precomputed.UW_rel_diff < 5e-8, ...
    'Precomputed big sparse with boundary postprocessing changed UW.');
assert(isfinite(result.boundary_precomputed.big_sparse_rel_res), ...
    'Boundary postprocessing with precomputed big sparse failed.');
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
opt.mob_big_sparse_chunk_pairs = 2;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.use_fmm = false;
end

function max_diff = compare_build_modes(q,opt)
[geom,basis] = build_test_data(q,opt);

geom_pre = geom;
geom_pre.opt.mob_big_sparse_build_mode = 'precomputed';
[big_pre,~] = buildMobPeanutBigSparseStokes(geom_pre,basis);

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
[big_stream,~] = buildMobPeanutBigSparseStokes(geom_stream,basis_stream);

diffs = zeros(4,1);
names = {'M_source_corr','M_pair_nonp','M_pair_proj','M_rbm_corr', ...
    'M_u_corr','M_u_cross','M_u_peanut'};
for k = 1:numel(names)
    diffs(k) = compare_optional_sparse(big_stream,big_pre,names{k});
end
max_diff = max(diffs);
end

function ram = compare_ram_estimates(q,opt)
n_pairs = count_close_pairs(q,opt.delta_pair,1);
N_check = ceil(opt.a_c*opt.N_c);

opt_pre = opt;
opt_pre.mob_big_sparse_build_mode = 'precomputed';
est_pre = estimateMobPeanutBigSparseRamStokes( ...
    numel(q),opt.N_c,N_check,n_pairs,opt_pre);

opt_pre_boundary = opt_pre;
opt_pre_boundary.get_bndry_field = 1;
est_pre_boundary = estimateMobPeanutBigSparseRamStokes( ...
    numel(q),opt.N_c,N_check,n_pairs,opt_pre_boundary);
assert(est_pre_boundary.retained_pair_basis_bytes > ...
        est_pre.retained_pair_basis_bytes, ...
    'Boundary precomputed estimate should include retained UB/YB pair factors.');

opt_stream = opt;
opt_stream.mob_big_sparse_build_mode = 'streaming';
est_stream = estimateMobPeanutBigSparseRamStokes( ...
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
F = [1 -0.5; -0.25 0.7; 0.4 -0.2; -1.15 0.0];
F = F - mean(F,1);
T = [0.2; -0.1; 0.05; -0.15];

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
F = [1 -0.5; -0.25 0.7; 0.4 -0.2; -1.15 0.0];
F = F - mean(F,1);
T = [0.2; -0.1; 0.05; -0.15];

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
