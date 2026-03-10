function [Q,lambda_proxy,it,gmres_tol,maxres] = solve_laplace_peanut_enhanced(q,g_body,delta_pair,N_peanut,visualise,gmres_tol,debug,use_fmm)
%SOLVE_LAPLACE_PEANUT_ENHANCED Solve exterior Dirichlet Laplace (a capacitance problem) with peanut-compressed 2B preconditioner.
%
% Syntax:
%   [Q,lambda_proxy,it,gmres_tol,maxres] = solve_laplace_peanut_enhanced(...)
%
% Inputs:
%   q          - Complex particle centers (P x 1).
%   g_body     - Constant boundary values per body (P x 1).
%   delta_pair - Pair threshold.
%   N_peanut   - Number of points on peanut boundary for pair compression.
%   visualise  - Plot diagnostics.
%   gmres_tol  - GMRES tolerance.
%   debug      - Build dense matrix diagnostics.
%   use_fmm    - Use fmm2d (of flatiron) for Laplace evaluations when available.
%
% Outputs:
%   Q           - Per-body unweighted sums of source strengths.
%   lambda_proxy- Compressed coarse source strengths.
%   it          - GMRES iteration count.
%   gmres_tol   - GMRES tolerance used.
%   maxres      - Max relative residual on independent boundary points.
%
% To test: call without inputs.
%
% Anna Broms, Mar 2026

if nargin==0
    test_solve_laplace_peanut;
    return
end

if nargin < 5 || isempty(visualise), visualise = 0; end
if nargin < 6 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 7 || isempty(debug), debug = false; end
if nargin < 8 || isempty(use_fmm), use_fmm = true; end

q = q(:);
g_body = g_body(:);
P = length(q);
assert(length(g_body)==P,'g_body must have one entry per particle.');

%% Parameters
maxit = 800;
solver_name = 'solve_laplace_peanut_enhanced';

opt = getLaplace2Dparams();
R = opt.alpha;

N_c = 80;
N_f = 150;
a_c = 1.2;
a_f = 1.2;

tol_c = 1e-10;
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = R*max([1-sep_c,0.01]);
Rp_f = R*max([1-sep_f,0.01]);

accstop = (R-Rp_c)^2/Rp_c;
if nargin < 3 || isempty(delta_pair)
    delta_pair = accstop;
end
if nargin < 4 || isempty(N_peanut)
    N_peanut = 400;
end

opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c;
opt.a_f = a_f;
opt.N_c = N_c;
opt.N_f = N_f;
opt.N_peanut = N_peanut;
opt.precomp = 1;
opt.pc = 1;
opt.delta_pair = delta_pair;
opt.P = P;
opt.Nclust = 100;
opt.cmap = 0;
opt.use_fmm = use_fmm;
opt.show_counter = true;
if visualise
    opt.visualse_grid = true; 
else
    opt.visualise_grid = false;
end
opt.rads = R*ones(P,1);

%% Build grids
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

rvec_in_c = zeros(P*N_c,1);
rout = zeros(P*nout,1);
for k = 1:P
    rvec_in_c((k-1)*N_c+1:k*N_c) = q(k)+rbase_in_c;
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
end

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

%% Basis factors
[UB_all,YB_all,UC_all,YC_all,Cmap] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt);
[UU,YY] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout]);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.opt = opt;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.rimage_vec = rimage_vec;
geom.rvec_in = rvec_in_c;

basis = struct();
basis.U = UU;
basis.Y = YY;
basis.Upf = UB_all;
basis.Ypf = YB_all;
basis.DC_all = UC_all;
basis.YC_all = YC_all;
basis.Cmap = Cmap;

%% RHS
fout = zeros(P*nout,1);
for k = 1:P
    fout((k-1)*nout+1:k*nout) = g_body(k);
end

%% Solve
if debug
    x = zeros(length(rout),1);
    CC = zeros(length(rout));
    for k = 1:length(rout)
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvec_laplace_peanut_enhanced(x,geom,basis);
    end
    figure(); imagesc(log10(abs(CC))); colorbar
    title([solver_name ': log_{10}|CC|'],'interpreter','none')
end

[tau,it,resvec,~] = helsing_gmres(@(x) matvec_laplace_peanut_enhanced(x,geom,basis), ...
    fout,length(rout),maxit,gmres_tol,1,rout);

figure(); semilogy(resvec)
title('GMRES convergence Laplace peanut enhanced','interpreter','latex')

%% Postprocess
n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+R*(cos(tb)+1i*sin(tb));
end

geom_eval = geom;
geom_eval.rcheck = rcheck_b;
[lam_c,lam_self,lam_f,lam_e,u_corr] = transform_laplace_peanut(tau,geom_eval,basis);

u_b = laplaceSingleLayerField(rvec_in_c,rcheck_b,lam_c,use_fmm) + u_corr;

g_true = zeros(P*n_bound,1);
for k = 1:P
    g_true((k-1)*n_bound+1:k*n_bound) = g_body(k);
end

maxres = max(abs(u_b-g_true))/max(1,max(abs(g_true)));
fprintf('Max surface relative residual at new nodes %.3e\n',maxres);

lambda_proxy = lam_c;

Q = zeros(P,1);
for k = 1:P
    idx = (k-1)*N_c+1:k*N_c;
    Q(k) = sum(lam_self(idx)) + sum(lam_f{k}) + sum(lam_e{k});
end

if visualise
    figure();
    plot(u_b); hold on; plot(g_true)
    title('Boundary values: lhs vs rhs (Laplace peanut)')

    figure();
    semilogy(abs(lambda_proxy))
    title('Compressed source strengths (Laplace peanut)')
end

end

function test_solve_laplace_peanut
fprintf('--- solve_laplace_peanut_enhanced self-test ---\n');

opt = getLaplace2Dparams();
R = opt.alpha;
% q = [0; 2*R+0.08*R; 6*R+1.5i*R];
% g_body = [1; -0.7; 0.25];

P = 40;
delta = 1e-2;
q = grow_cluster(P,delta,2,R);

g_body = rand(P,1);
delta_pair = 0.2;
N_peanut = 400;
visualise = 1; 

[Qp,~,itp,~,resp] = solve_laplace_peanut_enhanced(q,g_body,delta_pair,N_peanut,visualise,1e-10,0,true);
[Q2,~,it2,~,res2] = solve_laplace_2B_enhanced(q,g_body,delta_pair,0,1e-10,0,true);

fprintf('Peanut: it=%d, maxres=%.3e\n',itp,resp);
fprintf('2B    : it=%d, maxres=%.3e\n',it2,res2);
fprintf('Rel diff in Q (peanut vs 2B): %.3e\n',norm(Qp-Q2)/max(1,norm(Q2)));

end
