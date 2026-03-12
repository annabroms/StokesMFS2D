function [Q,lambda_all,it,gmres_tol,maxres] = solve_cap_1B(q,v_body,visualise,gmres_tol,debug,use_fmm,gmres_verbose)
%SOLVE_CAP_1B Solve exterior Dirichlet Laplace problem (capacitance: known
%voltages, unknown charges) with 1-body preconditioning.
%
% Syntax:
%   [Q,lambda_all,it,gmres_tol,maxres] = solve_cap_1B(...)
% Inputs:
%   q         - Complex particle centers (P x 1).
%   v_body    - Constant boundary values per body (P x 1).
%   visualise - Plot diagnostics.
%   gmres_tol - GMRES tolerance.
%   debug     - Build dense system matrix for diagnostics.
%   use_fmm   - Use fmm2d (of flatiron) for Laplace evaluations when available.
%
% Outputs:
%   Q          - Per-body sums of all source strengths belonging to each body.
%   lambda_all - Stacked source strengths on all one-body sources.
%   it         - GMRES iteration count.
%   gmres_tol  - GMRES tolerance used.
%   maxres     - Max relative residual on independent boundary points.
%
% Notes: what is solved is what is referred to as the modified exterior Laplace 
% BVP in the Stein & Barnett QFS paper from 2022. 
% The radius parameter is chosen with rad ~= 1 to avoid unit logarithmic
% capacity in 2D.
%
% To test: call without inputs.
%
% See also: solve_cap_2B, solve_cap_peanut, ...
%   solve_elast_1B, lapSLPfield.
%
% Anna Broms, Mar 2026

if nargin==0
    test_solve_cap_1B;
    return
end

if nargin < 3 || isempty(visualise), visualise = 0; end
if nargin < 4 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 5 || isempty(debug), debug = false; end
if nargin < 6 || isempty(use_fmm), use_fmm = true; end
if nargin < 7 || isempty(gmres_verbose), gmres_verbose = 0; end

q = q(:);
v_body = v_body(:);
P = numel(q);
assert(numel(v_body)==P,'v_body must have one entry per particle.');

maxit = 800;

if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

[geom,basis,~,R] = prepareLaplace1B(q,use_fmm);
opt = struct();
opt.gmres_verbose = gmres_verbose;

%% RHS
fout = zeros(length(geom.rout),1);
for k = 1:P
    fout(geom.target_ind{k}) = v_body(k);
end

%% Solve
if debug
    x = zeros(length(geom.rout),1);
    CC = zeros(length(geom.rout));
    ncols = length(geom.rout);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvec_laplace_1B(x,geom,basis);
    end
    figure(); imagesc(log10(abs(CC))); colorbar
    title([solver_name ': log_{10}|CC|'],'interpreter','none')
end

disp(' == Solving... == ');
[tau,it,resvec,~] = helsing_gmres(@(x) matvec_laplace_1B(x,geom,basis), ...
    fout,length(geom.rout),maxit,gmres_tol,opt,geom.rout);

figure(); semilogy(resvec)
title('GMRES convergence capacitance 1B','interpreter','latex')

disp(' == Postprocessing == ');
[lambda_all,lambda_body] = mapBoundaryToSources1B(tau,geom,basis);

%% Boundary check and extract net charges
n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+R*(cos(tb)+1i*sin(tb));
end

u_b = lapSLPfield(geom.rvec_in,rcheck_b,lambda_all,use_fmm);
b_true = zeros(P*n_bound,1);
for k = 1:P
    b_true((k-1)*n_bound+1:k*n_bound) = v_body(k);
end
maxres = max(abs(u_b-b_true))/max(1,max(abs(b_true)));
fprintf('Max surface relative residual at new nodes %.3e\n',maxres);

Q = zeros(P,1);
for k = 1:P
    Q(k) = sum(lambda_body{k});
end

if visualise
    figure();
    plot(u_b); hold on; plot(b_true)
    title('Boundary values: lhs vs rhs (Laplace 1B)')

    figure();
    semilogy(abs(lambda_all))
    title('Source strengths (proxy+enhancement) capacitance 1B')
    axis tight
end

end

function [lambda_all,lambda_body] = mapBoundaryToSources1B(tau,geom,basis)
P = length(geom.source_ind);
lambda_all = zeros(length(geom.rvec_in),1);
lambda_body = cell(P,1);

for k = 1:P
    tau_k = tau(geom.target_ind{k});
    lambda_k = basis.Y{k}*(basis.U{k}*tau_k);
    lambda_all(geom.source_ind{k}) = lambda_k;
    lambda_body{k} = lambda_k;
end
end

function res = matvec_laplace_1B(tau,geom,basis)
[lambda_all,lambda_body] = mapBoundaryToSources1B(tau,geom,basis);

res = lapSLPfield(geom.rvec_in,geom.rout,lambda_all,basis.use_fmm);

P = length(geom.target_ind);
for k = 1:P
    res(geom.target_ind{k}) = res(geom.target_ind{k}) - basis.Aii{k}*lambda_body{k};
end

res = res + tau;
end

function rcheck = buildExteriorPoints(q,R,npts)
q = q(:);
P = numel(q);
rcheck = zeros(0,1);

xmin = min(real(q)) - 3*R;
xmax = max(real(q)) + 3*R;
ymin = min(imag(q)) - 3*R;
ymax = max(imag(q)) + 3*R;

while numel(rcheck) < npts
    z = (xmin + (xmax-xmin)*rand(npts,1)) + 1i*(ymin + (ymax-ymin)*rand(npts,1));
    keep = true(size(z));
    for k = 1:P
        keep = keep & (abs(z-q(k)) > 1.2*R);
    end
    rcheck = [rcheck; z(keep)]; %#ok<AGROW>
    if numel(rcheck) > npts
        rcheck = rcheck(1:npts);
    end
end
end

function test_solve_cap_1B
fprintf('--- solve_cap_1B self-test ---\n');
close all; 
opt = getLaplace2Dparams();
R = opt.rad;
q = [0; 2*R+0.08*R; 6*R+1.5i*R];
P = numel(q);
v_body = [1; -0.7; 0.25];

[Q_it,lam_it,it_it,~,res_it] = solve_cap_1B(q,v_body,0,1e-10,0,true);
[geom,~,~,~] = prepareLaplace1B(q,false);

rhs = zeros(length(geom.rout),1);
for k = 1:P
    rhs(geom.target_ind{k}) = v_body(k);
end

A = lapSLPmat(geom.rvec_in,geom.rout);
[Yd,Ud] = getPseudoFactors(A,1e-14,0);
lam_dense = Yd*(Ud'*rhs);

Q_dense = zeros(P,1);
for k = 1:P
    Q_dense(k) = sum(lam_dense(geom.source_ind{k}));
end

n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+R*(cos(tb)+1i*sin(tb));
end
rcheck_ext = buildExteriorPoints(q,R,600);

u_it_b = lapSLPfield(geom.rvec_in,rcheck_b,lam_it,false);
u_dense_b = lapSLPfield(geom.rvec_in,rcheck_b,lam_dense,false);

u_it_ext = lapSLPfield(geom.rvec_in,rcheck_ext,lam_it,false);
u_dense_ext = lapSLPfield(geom.rvec_in,rcheck_ext,lam_dense,false);

rel_Q = norm(Q_it-Q_dense,inf)/max(1,norm(Q_dense,inf));
rel_b = norm(u_it_b-u_dense_b,inf)/max(1,norm(u_dense_b,inf));
rel_ext = norm(u_it_ext-u_dense_ext,inf)/max(1,norm(u_dense_ext,inf));

fprintf('Iterative: it=%d, maxres=%.3e\n',it_it,res_it);
fprintf('Dense-vs-iterative rel diff Q      : %.3e\n',rel_Q);
fprintf('Dense-vs-iterative rel diff bndry  : %.3e\n',rel_b);
fprintf('Dense-vs-iterative rel diff ext    : %.3e\n',rel_ext);

end
