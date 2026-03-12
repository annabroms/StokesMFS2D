function [v_body,lambda_all,it,gmres_tol,maxres] = solve_elast_2B(q,Q_body,delta_pair,visualise,gmres_tol,debug,use_fmm,gmres_verbose)
%SOLVE_ELAST_2B Solve exterior Laplace elastance problem (known charges, unknown voltages) with 2-body preconditioning.
%
% Syntax:
%   [v_body,lambda_all,it,gmres_tol,maxres] = solve_elast_2B(...)
%
% Inputs:
%   q          - Complex particle centers (P x 1).
%   Q_body     - Prescribed net charge per body (P x 1).
%   delta_pair - Pair threshold.
%   visualise  - Plot diagnostics.
%   gmres_tol  - GMRES tolerance.
%   debug      - Build dense system matrix for diagnostics.
%   use_fmm    - Use fmm2d (of flatiron) for Laplace evaluations when available.
%
% Outputs:
%   v_body    - Recovered constant boundary values per body (P x 1).
%   lambda_all- Stacked source strengths used in global field evaluation.
%   it        - GMRES iteration count.
%   gmres_tol - GMRES tolerance used.
%   maxres    - Max relative equipotential residual on independent boundary points.
%
% Notes:
%   The radius parameter is chosen with rad ~= 1 to avoid unit logarithmic
%   capacity in 2D.
%
% To test: call without inputs.
%
% See also: solve_elast_1B, solve_elast_peanut, solve_cap_2B, ...
%   getPairBasisLaplace, matvec_lap_2B_enhanced.
%
% Anna Broms, Mar 2026

if nargin==0
    test_solve_elast_2B;
    return
end

if nargin < 4 || isempty(visualise), visualise = 0; end
if nargin < 5 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 6 || isempty(debug), debug = false; end
if nargin < 7 || isempty(use_fmm), use_fmm = true; end
if nargin < 8 || isempty(gmres_verbose), gmres_verbose = 0; end

q = q(:);
Q_body = Q_body(:);
P = numel(q);
assert(numel(Q_body)==P,'Q_body must have one entry per particle.');

maxit = 800;

if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

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

opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c;
opt.a_f = a_f;
opt.N_c = N_c;
opt.N_f = N_f;
opt.N_peanut = 0;
opt.precomp = 1;
opt.pc = 1;
opt.delta_pair = delta_pair;
opt.P = P;
opt.Nclust = 100;
opt.use_fmm = use_fmm;
opt.show_counter = true;
opt.project_charge = true;
if visualise
    opt.visualise_grid = true;
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

rout = zeros(P*nout,1);
rvec_in_c = zeros(P*N_c,1);
coarse_source_ind = cell(P,1);
for k = 1:P
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
    idx = (k-1)*N_c+1:k*N_c;
    coarse_source_ind{k} = idx;
    rvec_in_c(idx) = q(k)+rbase_in_c;
end

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

%% Get 1- and 2-body basis 
[Upf,Ypf,~,~,~] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt);
[UU,YY] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout],true);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.opt = opt;
geom.rvec_out = rout;
geom.q = q;
geom.pairs = pairs;
geom.rimage_vec = rimage_vec;

basis = struct();
basis.U = UU;
basis.Y = YY;
basis.Upf = Upf;
basis.Ypf = Ypf;

[lambda0_c,u_rhs] = getChargeCompletionFlowLaplace(rvec_in_c,rout,coarse_source_ind,Q_body,use_fmm);

%% Solve
if debug
    x = zeros(length(rout),1);
    CC = zeros(length(rout));
    ncols = length(rout);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvec_lap_2B_enhanced(x,geom,basis,rout);
    end
    figure(); imagesc(log10(abs(CC))); colorbar
    title([solver_name ': log_{10}|CC|'],'interpreter','none')
end

disp(' == Solving... == ');
[tau,it,resvec,~] = helsing_gmres(@(x) matvec_lap_2B_enhanced(x,geom,basis,rout), ...
    u_rhs,length(rout),maxit,gmres_tol,opt,rout);

figure(); semilogy(resvec)
title('GMRES convergence elastance 2B','interpreter','latex')

disp(' == Postprocessing == ');
%% Postprocess
[rvec_in,coarse_ind,lambda_corr,~,~,~,~,lam_c_nonp,lam_f_nonp,lam_e_nonp] = ...
    getPairTransformationLaplace(tau,geom,basis);

lambda0_all = zeros(size(lambda_corr));
lambda0_all(1:P*N_c) = lambda0_c;
lambda_all = lambda0_all+lambda_corr;

n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+rad*(cos(tb)+1i*sin(tb));
end

u_b = lapSLPfield(rvec_in,rcheck_b,lambda_all,use_fmm);
v_body = zeros(P,1);
for k = 1:P
    idx = coarse_ind{k};
    lambda_tot_k = [-lam_c_nonp(idx); -lam_f_nonp{k}; -lam_e_nonp{k}];
    v_body(k) = sum(lambda_tot_k);
end
v_true = zeros(P*n_bound,1);
for k = 1:P
    v_true((k-1)*n_bound+1:k*n_bound) = v_body(k);
end
maxres = max(abs(u_b-v_true))/max(1,max(abs(v_true)));
fprintf('Max relative equipotential residual at new nodes %.3e\n',maxres);

if visualise
    figure();
    plot(u_b); hold on;
    plot(v_true);
    title('Boundary potential and per-body means (Laplace elastance 2B)')
end

end

function test_solve_elast_2B
fprintf('--- solve_elast_2B self-test ---\n');
close all;
run_two_way = false;

opt = getLaplace2Dparams();
rad = opt.rad;
delta = 1e-3;
visualise = 1; 
debug = 0; 
delta_pair = 0.2; 


rng(8);
P = 30;
q = grow_cluster(P,delta,2,rad);
Q_body = randn(P,1);
gmres_tol = 1e-7; 

[v2,~,it2,~,res2] = solve_elast_2B(q,Q_body,delta_pair,visualise,gmres_tol,debug,true);
[v1,~,it1,~,res1] = solve_elast_1B(q,Q_body,0,gmres_tol,0,true);

fprintf('2B: it=%d, maxres=%.3e\n',it2,res2);
fprintf('1B: it=%d, maxres=%.3e\n',it1,res1);
fprintf('Rel diff in v_body (2B vs 1B): %.3e\n',norm(v2-v1)/max(1,norm(v1)));

if run_two_way
    v_ref = randn(P,1);
    [Q_cap,~,~,~,~] = solve_cap_2B(q,v_ref,[],0,1e-10,0,true);
    [v_back,~,~,~,~] = solve_elast_2B(q,Q_cap,[],0,1e-10,0,true);
    rel_two = norm(v_back-v_ref,inf)/max(1,norm(v_ref,inf));
    fprintf('Two-way rel diff in v_body      : %.3e\n',rel_two);
end

alignfigs; 

end
