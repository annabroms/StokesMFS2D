clear;
close all;
set(0,'DefaultFigureVisible','off');

rng(13);

R = 2;
delta = 1e-2;
q_repeat = hexagonal_lattice(delta,1,R);
q_unique = [0:2*R+0.02:4*(2*R+0.02) (4*R+0.14)*1i]';

fprintf('=== Laplace Pair-Rotation Reuse Experiment (Mar 13, 2026) ===\n');
fprintf('Repeated-separation geometry: P=%d\n',numel(q_repeat));
run_pair_transform_checks(q_repeat,R);
run_solver_checks(q_repeat,R);

fprintf('\nUnique-separation geometry: P=%d\n',numel(q_unique));
run_group_count_check(q_unique,R);

function run_pair_transform_checks(q,R)
P = numel(q);
opt = getLaplace2Dparams(P,R);
opt.use_fmm = false;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.N_peanut = 200;
opt.delta_pair = 0.2;

configs = [ ...
    struct('precomp',1,'cmap',0), ...
    struct('precomp',0,'cmap',0), ...
    struct('precomp',1,'cmap',1)];

for kk = 1:numel(configs)
    cfg = configs(kk);
    opt.precomp = cfg.precomp;
    opt.cmap = cfg.cmap;

    opt_base = opt;
    opt_base.reuse_pair_basis_by_sep = false;
    opt_rot = opt;
    opt_rot.reuse_pair_basis_by_sep = true;

    [ctx_base,t_base] = build_context(q,opt_base);
    [ctx_rot,t_rot] = build_context(q,opt_rot);

    tau = randn(numel(ctx_base.geom2B.rvec_out),1);

    u_2b_base = matvec_lap_2B_enhanced(tau,ctx_base.geom2B,ctx_base.basis2B,ctx_base.geom2B.rvec_out);
    u_2b_rot = matvec_lap_2B_enhanced(tau,ctx_rot.geom2B,ctx_rot.basis2B,ctx_rot.geom2B.rvec_out);
    err_2b = relerr(u_2b_rot,u_2b_base);

    u_peanut_base = matvec_lap_peanut_enhanced(tau,ctx_base.geomPeanut,ctx_base.basisPeanut);
    u_peanut_rot = matvec_lap_peanut_enhanced(tau,ctx_rot.geomPeanut,ctx_rot.basisPeanut);
    err_peanut = relerr(u_peanut_rot,u_peanut_base);

    fprintf(['  precomp=%d cmap=%d: pairs=%d unique_seps=%d', ...
        ' setup %.3fs -> %.3fs, 2B matvec err %.3e, peanut matvec err %.3e\n'], ...
        cfg.precomp,cfg.cmap,size(ctx_base.geom2B.pairs,1),ctx_rot.basis2B.pair_cache.n_groups, ...
        t_base,t_rot,err_2b,err_peanut);
end
end

function run_solver_checks(q,R)
P = numel(q);
v_body = randn(P,1);
opt = getLaplace2Dparams(P,R);
opt.use_fmm = false;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.delta_pair = 0.2;
opt.gmres_tol = 1e-10;
opt.N_peanut = 100;

opt_2b = opt;
opt_2b.precomp = 1;

opt_2b_base = opt_2b;
opt_2b_base.reuse_pair_basis_by_sep = false;
opt_2b_rot = opt_2b;
opt_2b_rot.reuse_pair_basis_by_sep = true;

tic;
[Q_2b_base,sol_2b_base] = solve_cap_2B(q,v_body,opt_2b_base);
t_2b_base = toc;
tic;
[Q_2b_rot,sol_2b_rot] = solve_cap_2B(q,v_body,opt_2b_rot);
t_2b_rot = toc;

fprintf('  solve_cap_2B: relQ %.3e, relRes %.3e, time %.3fs -> %.3fs (x%.2f)\n', ...
    relerr(Q_2b_rot,Q_2b_base), ...
    relerr(sol_2b_rot.maxres,sol_2b_base.maxres), ...
    t_2b_base,t_2b_rot,t_2b_base/max(t_2b_rot,eps));

total_base = t_2b_base;
total_rot = t_2b_rot;

for cmap = 0:1
    opt_peanut = opt;
    opt_peanut.precomp = 1;
    opt_peanut.cmap = cmap;
    opt_peanut.compress_cmap = 0;

    opt_peanut_base = opt_peanut;
    opt_peanut_base.reuse_pair_basis_by_sep = false;
    opt_peanut_rot = opt_peanut;
    opt_peanut_rot.reuse_pair_basis_by_sep = true;

    tic;
    [Qp_base,solp_base] = solve_cap_peanut(q,v_body,opt_peanut_base);
    t_peanut_base = toc;
    tic;
    [Qp_rot,solp_rot] = solve_cap_peanut(q,v_body,opt_peanut_rot);
    t_peanut_rot = toc;

    fprintf('  solve_cap_peanut cmap=%d: relQ %.3e, relRes %.3e, time %.3fs -> %.3fs (x%.2f)\n', ...
        cmap,relerr(Qp_rot,Qp_base),relerr(solp_rot.maxres,solp_base.maxres), ...
        t_peanut_base,t_peanut_rot,t_peanut_base/max(t_peanut_rot,eps));

    total_base = total_base + t_peanut_base;
    total_rot = total_rot + t_peanut_rot;
end

fprintf('  combined solver time: %.3fs -> %.3fs (x%.2f)\n', ...
    total_base,total_rot,total_base/max(total_rot,eps));
end

function run_group_count_check(q,R)
P = numel(q);
opt = getLaplace2Dparams(P,R);
opt.use_fmm = false;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.reuse_pair_basis_by_sep = true;

[ctx,~] = build_context(q,opt);
fprintf('  pairs=%d unique_seps=%d\n',size(ctx.geom2B.pairs,1),ctx.basis2B.pair_cache.n_groups);
end

function [ctx,setup_time] = build_context(q,opt)
R = opt.rad;

P = numel(q);
N_c = opt.N_c;
N_f = opt.N_f;
a_c = opt.a_c;
a_f = opt.a_f;

tol_c = 1e-10;
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = getOptField(opt,'Rp_c',R*max([1-sep_c,0.01]));
Rp_f = getOptField(opt,'Rp_f',R*max([1-sep_f,0.01]));

nout = ceil(a_c*N_c);
tout = linspace(0,2*pi,nout+1)';
tout = tout(1:end-1);
rbase_out_c = R*(cos(tout)+1i*sin(tout));

tin = linspace(0,2*pi,N_c+1)';
tin = tin(1:end-1);
rbase_in_c = Rp_c*(cos(tin)+1i*sin(tin));

tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
rbase_in_f = Rp_f*(cos(tin_f)+1i*sin(tin_f));

tout_f = linspace(0,2*pi,ceil(a_f*N_f)+1)';
tout_f = tout_f(1:end-1);
rout_base_f = R*(cos(tout_f)+1i*sin(tout_f));

rvec_in_c = zeros(P*N_c,1);
rout = zeros(P*nout,1);
for k = 1:P
    idx = (k-1)*N_c+1:k*N_c;
    rvec_in_c(idx) = q(k)+rbase_in_c;
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
end

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

tic;
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_QV,pair_cache] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f,rimage_vec,refine,pairs,opt);
setup_time = toc;

[UU,YY] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout],false);

geom2B = struct();
geom2B.rbase_in_c = rbase_in_c;
geom2B.rbase_in_f = rbase_in_f;
geom2B.refine = refine;
geom2B.opt = opt;
geom2B.rvec_out = rout;
geom2B.q = q;
geom2B.pairs = pairs;
geom2B.rimage_vec = rimage_vec;
geom2B.pair_cache = pair_cache;

basis2B = struct();
basis2B.U = UU;
basis2B.Y = YY;
basis2B.Upf = UB_all;
basis2B.Ypf = YB_all;
basis2B.pair_cache = pair_cache;

geomPeanut = geom2B;
geomPeanut.rcheck = rout;
geomPeanut.rvec_in = rvec_in_c;

basisPeanut = basis2B;
basisPeanut.DC_all = UC_all;
basisPeanut.YC_all = YC_all;
basisPeanut.Cmap = Cmap;
basisPeanut.Cmap_QV = Cmap_QV;

ctx = struct();
ctx.geom2B = geom2B;
ctx.basis2B = basis2B;
ctx.geomPeanut = geomPeanut;
ctx.basisPeanut = basisPeanut;
end

function e = relerr(a,b)
e = norm(a-b,inf)/max(1,norm(b,inf));
end
