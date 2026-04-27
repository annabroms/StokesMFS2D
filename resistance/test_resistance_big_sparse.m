
% Compare loop and sparse resistance pair maps.

fprintf('--- test_resistance_big_sparse ---\n');

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

q = (2.06*(0:3)).';
opt = make_small_options(numel(q));

result = struct();
result.matvec_factored = compare_matvecs(q,opt,'factored');
result.matvec_combined = compare_matvecs(q,opt,'combined');
result.solve_factored = compare_solves(q,opt,'factored');
result.zero_pairs = compare_zero_close_pairs(opt);

fprintf('  PASSED\n');


function opt = make_small_options(P)
opt = get2Dparams(P,12,18);
opt.delta_pair = 0.12;
opt.N_peanut = 48;
opt.gmres_tol = 1e-7;
opt.maxit = 250;
opt.get_precomp_time = true;
opt.get_solve_time = true;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.use_fmm = false;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.res_big_sparse_chunk_pairs = 2;
opt.resistance = 1;
end

function result = compare_matvecs(q,opt,u_corr_mode)
opt.get_bndry_field = 1;
opt.res_big_sparse_u_corr_mode = u_corr_mode;
opt_loop = opt;
opt_loop.use_big_sparse = false;
[geom,basis] = build_test_data(q,opt_loop);

geom_loop = geom;
geom_loop.opt.use_big_sparse = false;

geom_sparse = geom;
geom_sparse.opt.use_big_sparse = true;
basis_sparse = basis;
[basis_sparse.big_sparse,stats] = buildResPeanutBigSparse( ...
    geom_sparse,basis_sparse);

rng(11);
tau = randn(2*numel(geom.rvec_out),1);
y_loop = matvec_res_peanut_enhanced(tau,geom_loop,basis);
y_sparse = matvec_peanut_big_sparse(tau,geom_sparse,basis_sparse);

result.rel_diff = relerr(y_sparse,y_loop);
result.big_sparse_stats = stats;
fprintf('  matvec u_mode=%s rel diff %.3e\n', ...
    u_corr_mode,result.rel_diff);
assert(result.rel_diff < 5e-10, ...
    'Sparse resistance matvec changed the solve-grid action.');
assert(strcmp(stats.u_corr_mode,u_corr_mode), ...
    'Sparse resistance build did not use the requested velocity-correction mode.');
if strcmp(u_corr_mode,'combined')
    assert(isfield(basis_sparse.big_sparse,'M_u_corr'), ...
        'Combined sparse resistance build did not create M_u_corr.');
else
    assert(isfield(basis_sparse.big_sparse,'M_u_cross') && ...
        isfield(basis_sparse.big_sparse,'M_u_peanut'), ...
        'Factored sparse resistance build did not create M_u_cross/M_u_peanut.');
end
end

function result = compare_solves(q,opt,u_corr_mode)
U = [1 -0.5; -0.25 0.7; 0.4 -0.2; -1.15 0.0];
W = [0.2; -0.1; 0.05; -0.15];

opt.res_big_sparse_u_corr_mode = u_corr_mode;
opt_loop = opt;
opt_loop.use_big_sparse = false;
opt_loop.get_bndry_field = 1;
[FT_loop,sol_loop] = solve_res_peanut_enhanced(q,U,W,opt_loop);

opt_sparse = opt;
opt_sparse.use_big_sparse = true;
opt_sparse.get_bndry_field = 1;
[FT_sparse,sol_sparse] = solve_res_peanut_enhanced(q,U,W,opt_sparse);

result.FT_rel_diff = relerr(FT_sparse,FT_loop);
assert(numel(sol_sparse.resvec) == numel(sol_loop.resvec), ...
    'Sparse resistance solve changed the residual-history length.');
result.resvec_diff = norm(sol_sparse.resvec-sol_loop.resvec,inf);
result.loop_boundary_rel_res = sol_loop.rel_res;
result.sparse_boundary_rel_res = sol_sparse.rel_res;
% The boundary residual itself is an accuracy diagnostic for this coarse
% test case; the regression check is that sparse and loop report the same
% diagnostic value.
result.rel_res_diff = abs(sol_sparse.rel_res-sol_loop.rel_res);
result.loop_iterations = sol_loop.it;
result.sparse_iterations = sol_sparse.it;

fprintf(['  solve u_mode=%s: it %d/%d, FT rel diff %.3e, ', ...
    'resvec diff %.3e, boundary rel res %.3e/%.3e, ', ...
    'boundary diff %.3e\n'], ...
    u_corr_mode,result.loop_iterations,result.sparse_iterations, ...
    result.FT_rel_diff,result.resvec_diff, ...
    result.loop_boundary_rel_res,result.sparse_boundary_rel_res, ...
    result.rel_res_diff);

assert(sol_sparse.it == sol_loop.it, ...
    'Sparse resistance solve changed the GMRES iteration count.');
assert(result.resvec_diff < 5e-10, ...
    'Sparse resistance solve changed the residual history.');
assert(result.FT_rel_diff < 5e-10, ...
    'Sparse resistance solve changed the force/torque output.');
assert(result.rel_res_diff < 5e-10, ...
    'Sparse resistance solve changed the boundary residual.');
assert(isfield(sol_sparse,'big_sparse_stats') && ...
    sol_sparse.big_sparse_stats.active, ...
    'Sparse resistance solve did not report active big_sparse_stats.');
assert(strcmp(sol_sparse.big_sparse_stats.u_corr_mode,u_corr_mode), ...
    'Sparse resistance solve did not report the requested velocity-correction mode.');
assert(~isfield(sol_sparse,'boundary_big_sparse_stats'), ...
    'Boundary postprocessing must not report a sparse boundary matrix.');
end

function result = compare_zero_close_pairs(opt)
q = [0; 4; 8];
opt0 = opt;
opt0.P = numel(q);
opt0.delta_pair = 0.1;
opt0.get_bndry_field = 0;

[geom,basis] = build_test_data(q,opt0);
geom_sparse = geom;
geom_sparse.opt.use_big_sparse = true;
[big_sparse,stats] = buildResPeanutBigSparse(geom_sparse,basis);
basis_sparse = basis;
basis_sparse.big_sparse = big_sparse;

rng(12);
tau = randn(2*numel(geom.rvec_out),1);
geom_loop = geom;
geom_loop.opt.use_big_sparse = false;
y_loop = matvec_res_peanut_enhanced(tau,geom_loop,basis);
y_sparse = matvec_peanut_big_sparse(tau,geom_sparse,basis_sparse);

result.n_pairs = stats.n_pairs;
result.rel_diff = relerr(y_sparse,y_loop);
fprintf('  zero-pair matvec rel diff %.3e\n',result.rel_diff);

assert(stats.n_pairs == 0, ...
    'Zero-close-pair test unexpectedly found close pairs.');
assert(strcmp(stats.reason,'no_close_pairs'), ...
    'Zero-close-pair sparse build should report no_close_pairs.');
assert(result.rel_diff < 5e-10, ...
    'Zero-close-pair sparse matvec changed the action.');
end

function [geom,basis] = build_test_data(q,opt)
q = q(:);
P = numel(q);
N_c = opt.N_c;
N_f = opt.N_f;
a_c = opt.a_c;

tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';
rbase_out_c = cos(tout_c)+1i*sin(tout_c);

tin_c = linspace(0,2*pi,N_c+1)';
tin_c = tin_c(1:end-1);
rbase_in_c = opt.Rp_c*cos(tin_c) + 1i*opt.Rp_c*sin(tin_c);

tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
rbase_in_f = opt.Rp_f*cos(tin_f) + 1i*opt.Rp_f*sin(tin_f);

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
[U,Y] = getSelfPseudo(1,rbase_in_c,rbase_out_c,[],[], ...
    [0,numel(rbase_out_c)],0,svd_opts);

opt_pair = opt;
opt_pair.project_force = false;
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU,pair_cache] = ...
    getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine, ...
    pairs,opt_pair,[],rbase_out_c,svd_opts);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.rimage_vec = rimage_vec;
geom.opt = opt_pair;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.rvec_in = rvec_in_c;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Upf = UB_all;
basis.Ypf = YB_all;
basis.DC_all = UC_all;
basis.YC_all = YC_all;
basis.Cmap = Cmap;
basis.Cmap_FU = Cmap_FU;
basis.pair_cache = pair_cache;
basis.Nii = stokSLPmat(rbase_in_c,rbase_out_c,1);
end

function e = relerr(a,b)
e = norm(a-b,inf)/max(1,norm(b,inf));
end
