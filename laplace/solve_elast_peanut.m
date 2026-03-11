function [v_body,lambda_proxy,it,gmres_tol,maxres] = solve_elast_peanut(q,Q_body,delta_pair,N_peanut,visualise,gmres_tol,debug,use_fmm,gmres_verbose)
%SOLVE_ELAST_PEANUT Solve exterior Laplace elastance problem (known charges, unknown voltages) with peanut-compressed 2B preconditioning.
%
% Syntax:
%   [v_body,lambda_proxy,it,gmres_tol,maxres] = solve_elast_peanut(...)
%
% Inputs:
%   q          - Complex particle centers (P x 1).
%   Q_body     - Prescribed net charge per body (P x 1).
%   delta_pair - Pair threshold.
%   N_peanut   - Number of points on peanut boundary for pair compression.
%   visualise  - Plot diagnostics.
%   gmres_tol  - GMRES tolerance.
%   debug      - Build dense system matrix diagnostics.
%   use_fmm    - Use fmm2d (of flatiron) for Laplace evaluations when available.
%
% Outputs:
%   v_body      - Recovered constant boundary values per body (P x 1).
%   lambda_proxy- Coarse proxy strengths after completion/correction merge.
%   it          - GMRES iteration count.
%   gmres_tol   - GMRES tolerance used.
%   maxres      - Max relative equipotential residual on independent boundary points.
%
% To test: call without inputs.
%
% See also: solve_elast_1B, solve_elast_2B, solve_cap_peanut, ...
%   transform_laplace_peanut, matvec_laplace_peanut_enhanced.
%
% Anna Broms, Mar 2026

if nargin==0
    test_solve_elast_peanut;
    return
end

if nargin < 5 || isempty(visualise), visualise = 0; end
if nargin < 6 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 7 || isempty(debug), debug = false; end
if nargin < 8 || isempty(use_fmm), use_fmm = true; end
if nargin < 9 || isempty(gmres_verbose), gmres_verbose = 0; end

q = q(:);
Q_body = Q_body(:);
P = numel(q);
assert(numel(Q_body)==P,'Q_body must have one entry per particle.');

maxit = 800;
solver_name = 'solve_elast_peanut';

opt = getLaplace2Dparams();
rad = opt.rad;
opt.gmres_verbose = gmres_verbose;

N_c = 80;
N_f = 150;
a_c = 1.2;
a_f = 1.2;

tol_c = 1e-10;
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = rad*max([1-sep_c,0.01]);
Rp_f = rad*max([1-sep_f,0.01]);

accstop = (rad-Rp_c)^2/Rp_c;
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
opt.project_charge = true;
if visualise
    opt.visualise_grid = true;
else
    opt.visualise_grid = false;
end
opt.rads = rad*ones(P,1);


%% Discretize
nout = ceil(a_c*N_c);
tout = linspace(0,2*pi,nout+1)';
tout = tout(1:end-1);
rbase_out_c = rad*(cos(tout)+1i*sin(tout));

tin = linspace(0,2*pi,N_c+1)';
tin = tin(1:end-1);
rbase_in_c = Rp_c*(cos(tin)+1i*sin(tin));

tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
rbase_in_f = Rp_f*(cos(tin_f)+1i*sin(tin_f));

rvec_in_c = zeros(P*N_c,1);
rout = zeros(P*nout,1);
coarse_source_ind = cell(P,1);
for k = 1:P
    idx = (k-1)*N_c+1:k*N_c;
    coarse_source_ind{k} = idx;
    rvec_in_c(idx) = q(k)+rbase_in_c;
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
end

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

%% Get 1- and 2-body basis functions
[UB_all,YB_all,UC_all,YC_all,Cmap] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt);
[UU,YY] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout],true);

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

%% Get rhs
[lambda0_c,u_rhs] = getChargeCompletionFlowLaplace(rvec_in_c,rout,coarse_source_ind,Q_body,use_fmm);

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
    u_rhs,length(rout),maxit,gmres_tol,opt,rout);

figure(); semilogy(resvec)
title('GMRES convergence elastance peanut','interpreter','latex')

%% Postprocess
n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+rad*(cos(tb)+1i*sin(tb));
end

geom_eval = geom;
geom_eval.rcheck = rcheck_b;
[lam_c,~,~,~,u_corr,~,lam_self_nonp,lam_f_nonp,lam_e_nonp] = ...
    transform_laplace_peanut(tau,geom_eval,basis);

u_b = laplaceSingleLayerField(rvec_in_c,rcheck_b,lambda0_c+lam_c,use_fmm);
u_b = u_b+u_corr;

v_body = zeros(P,1);
for k = 1:P
    idx = (k-1)*N_c+1:k*N_c;
    lambda_tot_k = [-lam_self_nonp(idx); -lam_f_nonp{k}; -lam_e_nonp{k}];
    v_body(k) = sum(lambda_tot_k);
end
v_true = zeros(P*n_bound,1);
for k = 1:P
    v_true((k-1)*n_bound+1:k*n_bound) = v_body(k);
end
maxres = max(abs(u_b-v_true))/max(1,max(abs(v_true)));
fprintf('Max relative equipotential residual at new nodes %.3e\n',maxres);

lambda_proxy = lambda0_c+lam_c;


if visualise
    figure();
    plot(u_b); hold on;
    plot(v_true);
    title('Boundary potential and per-body means (Laplace elastance peanut)')

    figure();
    semilogy(abs(lambda_proxy))
    title('Coarse proxy strengths (Laplace elastance peanut)')
end

end

function test_solve_elast_peanut
fprintf('--- solve_elast_peanut self-test ---\n');

close all; 
run_two_way = true;

opt = getLaplace2Dparams();
rad = opt.rad;
delta = 1e-3;
visualise = 1;
debug = 0; % determine system matrix by using matvec with columns of the identity matrix as input

rng(9);
P = 20;
q = grow_cluster(P,delta,2,rad);
Q_body = randn(P,1);

delta_pair = 0.2;
N_peanut = 400;

[vp,~,itp,~,resp] = solve_elast_peanut(q,Q_body,delta_pair,N_peanut,visualise,1e-10,debug,true);
[v2,~,it2,~,res2] = solve_elast_2B(q,Q_body,delta_pair,0,1e-10,0,true);

fprintf('Peanut: it=%d, maxres=%.3e\n',itp,resp);
fprintf('2B    : it=%d, maxres=%.3e\n',it2,res2);
fprintf('Rel diff in v_body (peanut vs 2B): %.3e\n',norm(vp-v2)/max(1,norm(v2)));

if run_two_way
    if visualise
        disp('Press any key to continue')
        pause();
    end
    v_ref = randn(P,1);
    [Q_cap,~,~,~,~] = solve_cap_peanut(q,v_ref,delta_pair,N_peanut,0,1e-10,0,true);
    [v_back,~,~,~,~] = solve_elast_peanut(q,Q_cap,delta_pair,N_peanut,0,1e-10,0,true);
    rel_two = norm(v_back-v_ref,inf)/max(1,norm(v_ref,inf));
    fprintf('Two-way rel diff in v_body         : %.3e\n',rel_two);
end

end
