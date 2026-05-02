function result = test_mobility_sparse_map_coarse
%TEST_MOBILITY_SPARSE_MAP_COARSE Small parity checks for serial vs parallel pair precompute.

fprintf('--- test_mobility_sparse_map_coarse ---\n');

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

P = 10;
phi = 0.65;
rng_seed = 240424;
[q,meta] = build_small_geometry(P,phi,rng_seed);
opt = make_small_options(P);

fprintf('  geometry: P=%d, phi=%.3f, close pairs=%d\n', ...
    P,meta.phi,count_close_pairs(q,opt.delta_pair,1));

result = struct();
result.cases = repmat(struct(),0,1);
case_id = 0;
for direct_u_corr = [true false]
    for sparse_map_coarse = [false true]
        case_id = case_id + 1;
        result.cases(case_id) = compare_case( ...
            q,opt,direct_u_corr,sparse_map_coarse);
    end
end

assert(max([result.cases.transform_rel_diff]) < 5e-10, ...
    'Parallel precompute changed the mobility sparse transform.');
assert(max([result.cases.matvec_rel_diff]) < 5e-10, ...
    'Parallel precompute changed the mobility sparse matvec.');
assert(max([result.cases.solve_rel_diff]) < 5e-8, ...
    'Parallel precompute changed the mobility solve output.');
assert(max([result.cases.resvec_overlap_diff]) < 5e-10, ...
    'Parallel precompute changed the mobility GMRES residual history.');
assert(max([result.cases.final_resvec_diff]) < 5e-10, ...
    'Parallel precompute changed the final mobility GMRES residual.');
assert(all([result.cases.it_serial] == [result.cases.it_parallel]), ...
    'Parallel precompute changed the mobility GMRES iteration count.');

fprintf('  PASSED\n');
end

function case_result = compare_case(q,opt,direct_u_corr,sparse_map_coarse)
ensure_pool(2);

opt_serial = opt;
opt_serial.use_big_sparse = true;
opt_serial.mob_sparse_map_coarse = sparse_map_coarse;
opt_serial.big_sparse_direct_u_corr = direct_u_corr;
opt_serial.parallel_precomp = false;

opt_parallel = opt_serial;
opt_parallel.parallel_precomp = true;
opt_parallel.parallel_precomp_chunk_pairs = 2;

[geom_serial,basis_serial] = build_test_data(q,opt_serial);
[geom_parallel,basis_parallel] = build_test_data(q,opt_parallel);

[basis_serial.big_sparse,stats_serial] = buildMobPeanutBigSparse( ...
    geom_serial,basis_serial);
[basis_parallel.big_sparse,stats_parallel] = buildMobPeanutBigSparse( ...
    geom_parallel,basis_parallel);

rng(11);
tau = randn(2*numel(geom_serial.rvec_out),1);
[lam_serial_x,lam_serial_y,lambda_self_serial,u_corr_serial] = ...
    transform_mob_peanut_big_sparse_stokes(tau,geom_serial,basis_serial);
[lam_parallel_x,lam_parallel_y,lambda_self_parallel,u_corr_parallel] = ...
    transform_mob_peanut_big_sparse_stokes(tau,geom_parallel,basis_parallel);

y_serial = matvec_peanut_big_sparse(tau,geom_serial,basis_serial);
y_parallel = matvec_peanut_big_sparse(tau,geom_parallel,basis_parallel);

F = randn(numel(q),2);
F = F - mean(F,1);
T = randn(numel(q),1);
[UW_serial,sol_serial] = solve_mob_peanut_enhanced(q,F,T,opt_serial);
[UW_parallel,sol_parallel] = solve_mob_peanut_enhanced(q,F,T,opt_parallel);

case_result = struct();
case_result.direct_u_corr = direct_u_corr;
case_result.sparse_map_coarse = sparse_map_coarse;
case_result.transform_rel_diff = max([ ...
    relerr(lam_parallel_x,lam_serial_x), ...
    relerr(lam_parallel_y,lam_serial_y), ...
    relerr(lambda_self_parallel,lambda_self_serial), ...
    relerr(u_corr_parallel,u_corr_serial)]);
case_result.matvec_rel_diff = relerr(y_parallel,y_serial);
case_result.solve_rel_diff = relerr(UW_parallel,UW_serial);
case_result.resvec_overlap_diff = compare_vec_overlap( ...
    sol_parallel.resvec,sol_serial.resvec);
case_result.final_resvec_diff = compare_final_entry( ...
    sol_parallel.resvec,sol_serial.resvec);
case_result.it_serial = sol_serial.it;
case_result.it_parallel = sol_parallel.it;

fprintf(['  direct_u_corr=%d mob_sparse_map_coarse=%d transform %.3e ', ...
    'matvec %.3e solve %.3e it %d/%d resvec %.3e final %.3e\n'], ...
    direct_u_corr,sparse_map_coarse,case_result.transform_rel_diff, ...
    case_result.matvec_rel_diff,case_result.solve_rel_diff, ...
    case_result.it_serial,case_result.it_parallel, ...
    case_result.resvec_overlap_diff,case_result.final_resvec_diff);
end

function [q,meta] = build_small_geometry(P,phi,rng_seed)
geom_opt = struct('domain','boxed','phi',phi,'rad',1, ...
    'min_gap',1e-3,'n_sweeps',200,'rng_seed',rng_seed, ...
    'visualise',false);
[q,meta] = random_discs_mc(P,geom_opt);
end

function opt = make_small_options(P)
opt = get2Dparams(P,12,18);
opt.a_c = 1.25;
opt.delta_pair = 0.2;
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
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.parallel_precomp_chunk_pairs = 2;
opt.use_fmm = false;
opt.mob_big_sparse_chunk_pairs = 2;
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

svd_opts = struct( ...
    'column_weight',logical(getOptField(opt,'column_weight',false)), ...
    'left_weight',logical(getOptField(opt,'left_weight',false)));
[U,Y,Lc] = getSelfPseudoMobilityStokes( ...
    1,q,rbase_in_c,rbase_out_c,[],[0,ceil(a_c*N_c)],svd_opts);

opt_pair = opt;
opt_pair.project_force = true;
opt_pair.project = true;
opt_pair.pair_basis_debug = false;
opt_pair.rad = ones(P,1);
if logical(opt_pair.parallel_precomp)
    opt_pair.parallel_precomp_chunk_pairs = 2;
end
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU,pair_cache] = ...
    getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine, ...
    pairs,opt_pair,Lc{1},rbase_out_c,svd_opts);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.rvec_in = rvec_in_c;
geom.rimage_vec = rimage_vec;
geom.opt = opt_pair;
geom.opt.get_bndry_field = 0;
geom.opt.parallel_solve = false;
geom.opt.resistance = 0;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.refine = refine;

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

function e = relerr(a,b)
den = max(1,norm(b(:),inf));
e = norm(a(:)-b(:),inf)/den;
end

function value = compare_vec_overlap(a,b)
a = a(:);
b = b(:);
n = min(numel(a),numel(b));
if n == 0
    value = 0;
    return
end
value = norm(a(1:n)-b(1:n),inf);
end

function value = compare_final_entry(a,b)
a = a(:);
b = b(:);
assert(~isempty(a) && ~isempty(b), ...
    'Residual history is empty in mobility sparse-map comparison.');
value = abs(a(end) - b(end));
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
