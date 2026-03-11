function [Q,lambda_all,it,gmres_tol,maxres] = solve_cap_2B(q,v_body,delta_pair,visualise,gmres_tol,debug,use_fmm)
%SOLVE_CAP_2B Solve exterior Dirichlet Laplace problem
%(capacitance: known voltages, unknown charges) with 2 body preconditioning
%(no compression).
%
% Syntax:
%   [Q,lambda_all,it,gmres_tol,maxres] = solve_cap_2B(...)
% Inputs:
%   q          - Complex particle centers (P x 1).
%   v_body     - Constant boundary values per body (P x 1).
%   delta_pair - Pair threshold.
%   visualise  - Plot diagnostics.
%   gmres_tol  - GMRES tolerance.
%   debug      - Build dense system matrix for diagnostics.
%   use_fmm    - Use fmm2d (of flatiron) for Laplace evaluations when available.
%
% Outputs:
%   Q          - Per-body unweighted sums of source strengths.
%   lambda_all - Stacked source strengths used in global field evaluation.
%   it         - GMRES iteration count.
%   gmres_tol  - GMRES tolerance used.
%   maxres     - Max relative residual on independent boundary points.
%
% To test: call without inputs.
%
% See also: solve_cap_1B, solve_cap_peanut, solve_elast_2B, ...
%   getPairBasisLaplace, matvec_lap_2B_enhanced.
%
% Anna Broms, Mar 2026

if nargin==0
    test_solve_cap_2B;
    return
end

if nargin < 4 || isempty(visualise), visualise = 0; end
if nargin < 5 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 6 || isempty(debug), debug = false; end
if nargin < 7 || isempty(use_fmm), use_fmm = true; end

q = q(:);
v_body = v_body(:);
P = numel(q);

assert(numel(v_body)==P,'v_body must have one entry per particle.');

%% Parameters
maxit = 800;
solver_name = 'cap_2B';

opt = getLaplace2Dparams();
R = opt.rad;

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

if visualise
    opt.visualise_grid = true;
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

rout = zeros(P*nout,1);
for k = 1:P
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
end

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

%% Basis factors
[Upf,Ypf,~,~,~] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt);
[UU,YY] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout]);

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

%% RHS
fout = zeros(P*nout,1);
for k = 1:P
    fout((k-1)*nout+1:k*nout) = v_body(k);
end

%% Solve
if debug
    x = zeros(length(rout),1);
    CC = zeros(length(rout));
    for k = 1:length(rout)
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvec_lap_2B_enhanced(x,geom,basis,rout);
    end
    figure(); imagesc(log10(abs(CC))); colorbar
    title([solver_name ': log_{10}|CC|'],'interpreter','none')
end

[tau,it,resvec,~] = helsing_gmres(@(x) matvec_lap_2B_enhanced(x,geom,basis,rout), ...
    fout,length(rout),maxit,gmres_tol,1,rout);

figure(); semilogy(resvec)
title('GMRES convergence Laplace 2B enhanced','interpreter','latex')

%% Postprocess
[rvec_in,coarse_ind,lambda_all,lam_c,lam_f,lam_e] = getPairTransformationLaplace(tau,geom,basis);

n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+R*(cos(tb)+1i*sin(tb));
end

u_b = laplaceSingleLayerField(rvec_in,rcheck_b,lambda_all,use_fmm);
g_true = zeros(P*n_bound,1);
for k = 1:P
    g_true((k-1)*n_bound+1:k*n_bound) = v_body(k);
end

maxres = max(abs(u_b-g_true))/max(1,max(abs(g_true)));
fprintf('Max surface relative residual at new nodes %.3e\n',maxres);

% Determine net charges for all objects
Q = zeros(P,1);
for k = 1:P
    idx = coarse_ind{k};
    Q(k) = sum(lam_c(idx)) + sum(lam_f{k}) + sum(lam_e{k});
end

if visualise
    figure();
    plot(u_b); hold on; plot(g_true)
    title('Boundary values: lhs vs rhs (Laplace 2B)')
end

end

function test_solve_cap_2B
fprintf('--- solve_cap_2B self-test ---\n');

close all; 
opt = getLaplace2Dparams();
R = opt.rad;
delta = 0.001; 
% q = [0; 2*R+delta*R; 6*R+1.5i*R];
% v_body = [1; -0.7; 0.25];

rng(8);
P = 40;
q = grow_cluster(P,delta,2,R);
v_body = rand(P,1); 
visualise = 1; 

[Q2,~,it2,~,res2] = solve_cap_2B(q,v_body,[],visualise,1e-10,0,true);
[Q1,~,it1,~,res1] = solve_cap_1B(q,v_body,visualise,1e-10,0,true);

fprintf('2B: it=%d, maxres=%.3e\n',it2,res2);
fprintf('1B: it=%d, maxres=%.3e\n',it1,res1);
fprintf('Rel diff in Q (2B vs 1B): %.3e\n',norm(Q2-Q1)/max(1,norm(Q1)));

end
