function [v_body,lambda_all,it,gmres_tol,maxres] = solve_elast_1B(q,Q_body,visualise,gmres_tol,debug,use_fmm,gmres_verbose)
%SOLVE_ELAST_1B Solve exterior Laplace elastance problem (known charges, unknown voltages) with 1-body preconditioning.
%
% Syntax:
%   [v_body,lambda_all,it,gmres_tol,maxres] = solve_elast_1B(...)
%
% Inputs:
%   q         - Complex particle centers (P x 1).
%   Q_body    - Prescribed net charge per body (P x 1).
%   visualise - Plot diagnostics.
%   gmres_tol - GMRES tolerance.
%   debug     - Build dense system matrix for diagnostics.
%   use_fmm   - Use fmm2d (of flatiron) for Laplace evaluations when available.
%
% Outputs:
%   v_body    - Recovered constant boundary values per body (P x 1).
%   lambda_all- Stacked source strengths (all one-body sources).
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
% See also: solve_elast_2B, solve_elast_peanut, solve_cap_1B, ...
%   getChargeCompletionFlowLaplace, lapSLPField.
%
% Anna Broms, Mar 2026

if nargin==0
    test_solve_elast_1B;
    return
end

if nargin < 3 || isempty(visualise), visualise = 0; end
if nargin < 4 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 5 || isempty(debug), debug = false; end
if nargin < 6 || isempty(use_fmm), use_fmm = true; end
if nargin < 7 || isempty(gmres_verbose), gmres_verbose = 0; end

q = q(:);
Q_body = Q_body(:);
P = numel(q);
assert(numel(Q_body)==P,'Q_body must have one entry per particle.');

maxit = 800;

if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

%% Build geometry and basis functions
[geom,basis,~,rad] = prepareLaplace1B(q,use_fmm,visualise,true);
opt = struct();
opt.gmres_verbose = gmres_verbose;

%% Build rhs
[lambda0,u_rhs] = getChargeCompletionFlowLaplace(geom.rvec_in,geom.rout,geom.source_ind,Q_body,use_fmm);

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
        CC(:,k) = matvec_elast_1B(x,geom,basis);
    end
    figure(); imagesc(log10(abs(CC))); colorbar
    title([solver_name ': log_{10}|CC|'],'interpreter','none')
end

disp(' == Solving... == ');
[tau,it,resvec,~] = helsing_gmres(@(x) matvec_elast_1B(x,geom,basis), ...
    u_rhs,length(geom.rout),maxit,gmres_tol,opt,geom.rout);

figure(); semilogy(resvec)
title('GMRES convergence elastance 1B','interpreter','latex')

disp(' == Postprocessing == ');
%% Postprocess

[lambda_compute,lambda_body] = mapBoundaryToSourcesElast1B(tau,geom,basis);
lambda_all = lambda0+lambda_compute; 
% these are the total point charges. Have to use the projected sources here

% Create check surface on everyone for evaluation of boundary voltages
n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+rad*(cos(tb)+1i*sin(tb));
end

%Evaluate lhs, based on the representation 
u_b = lapSLPField(geom.rvec_in,rcheck_b,lambda_all,use_fmm);

% Use unprojected sources to determine net voltages
v_body = zeros(P,1);
for k = 1:P
    v_body(k) = -sum(lambda_body{k});
end

% the net charges go into the rhs, that should equal the lhs
v_bnd = zeros(P*n_bound,1);
for k = 1:P
    v_bnd((k-1)*n_bound+1:k*n_bound) = v_body(k);
end
maxres = max(abs(u_b-v_bnd))/max(1,max(abs(v_bnd)));
fprintf('Max relative equipotential residual at new nodes %.3e\n',maxres);


if visualise
    figure();
    plot(u_b); hold on;
    plot(v_bnd);
    title('Boundary potential and per-body means (Laplace elastance 1B)')
end

end

function [lambda_all,lambda_body] = mapBoundaryToSourcesElast1B(tau,geom,basis)
P = length(geom.source_ind);
lambda_all = zeros(length(geom.rvec_in),1);
lambda_body = cell(P,1);

for k = 1:P
    tau_k = tau(geom.target_ind{k});
    lambda_k = basis.Y{k}*(basis.U{k}*tau_k);
    lambda_body{k} = lambda_k; % store unprojected sources, needed to evaluate voltages
    lambda_k = lambda_k - mean(lambda_k); %project off charge contribution
    lambda_all(geom.source_ind{k}) = lambda_k;
    
end
end

function res = matvec_elast_1B(tau,geom,basis)
[lambda_all,lambda_body] = mapBoundaryToSourcesElast1B(tau,geom,basis);

res = lapSLPField(geom.rvec_in,geom.rout,lambda_all,basis.use_fmm);

P = length(geom.target_ind);
for k = 1:P
    res(geom.target_ind{k}) = res(geom.target_ind{k})-basis.Aii{k}*lambda_body{k};
end

res = res+tau;
end

function rcheck = buildExteriorPoints(q,rad,npts)
q = q(:);
P = numel(q);
rcheck = zeros(0,1);

xmin = min(real(q))-3*rad;
xmax = max(real(q))+3*rad;
ymin = min(imag(q))-3*rad;
ymax = max(imag(q))+3*rad;

while numel(rcheck) < npts
    z = (xmin + (xmax-xmin)*rand(npts,1)) + 1i*(ymin + (ymax-ymin)*rand(npts,1));
    keep = true(size(z));
    for k = 1:P
        keep = keep & (abs(z-q(k)) > 1.2*rad);
    end
    rcheck = [rcheck; z(keep)]; %#ok<AGROW>
    if numel(rcheck) > npts
        rcheck = rcheck(1:npts);
    end
end
end

function test_solve_elast_1B
fprintf('--- solve_elast_1B self-test ---\n');
close all;
run_two_way = true;
visualise = false; 
debug = false; % track system matrix and visualise

opt = getLaplace2Dparams();
rad = opt.rad;
delta = 1e-3; 
q = [0; rad*(2+delta)];
P = numel(q);
Q_body = [1; -2];

[v_it,lam_it,it_it,~,res_it] = solve_elast_1B(q,Q_body,visualise,1e-10,debug,true);
[geom,~,~,~] = prepareLaplace1B(q,false,visualise,true);
[lambda0,u_rhs] = getChargeCompletionFlowLaplace(geom.rvec_in,geom.rout,geom.source_ind,Q_body,false);


A_dense = buildDenseElastSystemMatrix(geom);
[Yd,Ud] = getPseudoFactors(A_dense,1e-14,0);
lambda_unproj_dense = Yd*(Ud'*u_rhs);
lambda_corr_dense = projectChargePerBody(lambda_unproj_dense,geom.source_ind);
lam_dense = lambda0+lambda_corr_dense;

n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+rad*(cos(tb)+1i*sin(tb));
end
rcheck_ext = buildExteriorPoints(q,rad,600);

u_it_b = lapSLPField(geom.rvec_in,rcheck_b,lam_it,false);
u_dense_b = lapSLPField(geom.rvec_in,rcheck_b,lam_dense,false);
u_it_ext = lapSLPField(geom.rvec_in,rcheck_ext,lam_it,false);
u_dense_ext = lapSLPField(geom.rvec_in,rcheck_ext,lam_dense,false);

v_dense = zeros(P,1);
for k = 1:P
    idx = geom.source_ind{k};
    v_dense(k) = -sum(lambda_unproj_dense(idx));
end

rel_v = norm(v_it-v_dense,inf)/max(1,norm(v_dense,inf));
rel_b = norm(u_it_b-u_dense_b,inf)/max(1,norm(u_dense_b,inf));
rel_ext = norm(u_it_ext-u_dense_ext,inf)/max(1,norm(u_dense_ext,inf));

fprintf('1-body iterative: it=%d, maxres=%.3e\n',it_it,res_it);
fprintf('Dense-vs-iterative rel diff v_body : %.3e\n',rel_v);
fprintf('Dense-vs-iterative rel diff bndry  : %.3e\n',rel_b);
fprintf('Dense-vs-iterative rel diff ext    : %.3e\n',rel_ext);

if run_two_way
    disp('Press key to continue...')
    pause();
    v_ref = rand(P,1); 
    [Q_cap,~,~,~,~] = solve_cap_1B(q,v_ref,0,1e-10,debug,true);
    [v_back,~,~,~,~] = solve_elast_1B(q,Q_cap,0,1e-10,0,true);
    rel_two = norm(v_back-v_ref,inf)/max(1,norm(v_ref,inf));
    fprintf('Two-way rel diff in v_body         : %.3e\n',rel_two);
end

function A_dense = buildDenseElastSystemMatrix(geom)
S = lapSLPmat(geom.rvec_in,geom.rout);
A_dense = S;

P = length(geom.source_ind);
for k = 1:P
    cols = geom.source_ind{k};
    nk = length(cols);
    if nk==0
        continue
    end
    row_mean = sum(A_dense(:,cols),2)/nk;
    A_dense(:,cols) = A_dense(:,cols)-row_mean*ones(1,nk);
end

for k = 1:P
    rows = geom.target_ind{k};
    cols = geom.source_ind{k};
    if isempty(rows) || isempty(cols)
        continue
    end
    A_dense(rows,cols) = A_dense(rows,cols)+1;
end
end

function lambda_proj = projectChargePerBody(lambda_unproj,source_ind)
lambda_proj = lambda_unproj;
P = length(source_ind);
for k = 1:P
    idx = source_ind{k};
    if isempty(idx)
        continue
    end
    lambda_proj(idx) = lambda_proj(idx)-mean(lambda_proj(idx));
end
end

end
