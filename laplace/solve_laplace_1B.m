function [Q,lambda_all,it,gmres_tol,maxres] = solve_laplace_1B(q,g_body,visualise,gmres_tol,debug,use_fmm)
%SOLVE_LAPLACE_1B Solve exterior Dirichlet Laplace problem (capacitance) with 1-body preconditioning.
%
% Syntax:
%   [Q,lambda_all,it,gmres_tol,maxres] = solve_laplace_1B(...)
%
% Inputs:
%   q         - Complex particle centers (P x 1).
%   g_body    - Constant boundary values per body (P x 1).
%   visualise - Plot diagnostics.
%   gmres_tol - GMRES tolerance.
%   debug     - Build dense system matrix for diagnostics.
%   use_fmm   - Use fmm2d (flatiron) for Laplace evaluations when available.
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
%
% To test: call without inputs.
%
% Anna Broms, Mar 2026

if nargin==0
    test_solve_laplace_1B;
    return
end

if nargin < 3 || isempty(visualise), visualise = 0; end
if nargin < 4 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 5 || isempty(debug), debug = false; end
if nargin < 6 || isempty(use_fmm), use_fmm = true; end

q = q(:);
g_body = g_body(:);
P = length(q);
assert(length(g_body)==P,'g_body must have one entry per particle.');

maxit = 800;
solver_name = 'solve_laplace_1B';

[geom,basis,~,R] = prepareLaplace1B(q,use_fmm);

%% RHS
fout = zeros(length(geom.rout),1);
for k = 1:P
    fout(geom.target_ind{k}) = g_body(k);
end

%% Solve
if debug
    x = zeros(length(geom.rout),1);
    CC = zeros(length(geom.rout));
    for k = 1:length(geom.rout)
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvec_laplace_1B(x,geom,basis);
    end
    figure(); imagesc(log10(abs(CC))); colorbar
    title([solver_name ': log_{10}|CC|'],'interpreter','none')
end

[tau,it,resvec,~] = helsing_gmres(@(x) matvec_laplace_1B(x,geom,basis), ...
    fout,length(geom.rout),maxit,gmres_tol,1,geom.rout);

figure(); semilogy(resvec)
title('GMRES convergence Laplace 1B','interpreter','latex')

[lambda_all,lambda_body] = mapBoundaryToSources1B(tau,geom,basis);

%% Boundary check and extract net charges
n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+R*(cos(tb)+1i*sin(tb));
end

u_b = laplaceSingleLayerField(geom.rvec_in,rcheck_b,lambda_all,use_fmm);
g_true = zeros(P*n_bound,1);
for k = 1:P
    g_true((k-1)*n_bound+1:k*n_bound) = g_body(k);
end
maxres = max(abs(u_b-g_true))/max(1,max(abs(g_true)));
fprintf('Max surface relative residual at new nodes %.3e\n',maxres);

Q = zeros(P,1);
for k = 1:P
    Q(k) = sum(lambda_body{k});
end

if visualise
    figure();
    plot(u_b); hold on; plot(g_true)
    title('Boundary values: lhs vs rhs (Laplace 1B)')

    figure();
    semilogy(abs(lambda_all))
    title('Source strengths (Laplace 1B)')
end

end

function [geom,basis,opt,R] = prepareLaplace1B(q,use_fmm)
P = length(q);

opt = getLaplace2Dparams();
R = opt.alpha;

N_c = 80;
N_f = 150;
a_c = 1.2;
a_f = 1.2;

% Keep same proxy-separation rule, but with scaled radii.
tol_c = 1e-10;
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = R*max([1-sep_c,0.01]);
Rp_f = R*max([1-sep_f,0.01]);

accstop = (R-Rp_c)^2/Rp_c;

opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c;
opt.a_f = a_f;
opt.N_c = N_c;
opt.N_f = N_f;
opt.N_peanut = 0;
opt.precomp = 1;
opt.pc = 0;
opt.delta_pair = accstop;
opt.P = P;
opt.Nclust = 100;
opt.cmap = 0;
opt.use_fmm = use_fmm;
opt.show_counter = false;
opt.visualise_grid = false;
opt.rads = R*ones(P,1);

nout = ceil(a_c*N_c);
tout = linspace(0,2*pi,nout+1)';
tout = tout(1:end-1);
rbase_out_c = R*(cos(tout)+1i*sin(tout));

tin = linspace(0,2*pi,N_c+1)';
tin = tin(1:end-1);
rbase_in_c = Rp_c*(cos(tin)+1i*sin(tin));

[cent_clust_cells,~,coll_clust_cells,~,~,pairs] = getEnhancedGrid(q,opt);

rin_body = cell(P,1);
rout_body = cell(P,1);
source_ind = cell(P,1);
target_ind = cell(P,1);
Aii = cell(P,1);
U = cell(P,1);
Y = cell(P,1);

r_chunks = cell(P,1);
t_chunks = cell(P,1);
source_start = 1;
target_start = 1;

for k = 1:P
    rin_k = [q(k)+rbase_in_c; cent_clust_cells{k}];
    rout_k = [q(k)+rbase_out_c; coll_clust_cells{k}];

    rin_body{k} = rin_k;
    rout_body{k} = rout_k;

    ns = length(rin_k);
    nt = length(rout_k);

    source_ind{k} = source_start:source_start+ns-1;
    target_ind{k} = target_start:target_start+nt-1;
    source_start = source_start+ns;
    target_start = target_start+nt;

    Aii{k} = lapSLPmat(rin_k,rout_k);
    [Yk,Uk] = getPseudoFactors(Aii{k},1e-14,0);
    U{k} = Uk';
    Y{k} = Yk;

    r_chunks{k} = rin_k;
    t_chunks{k} = rout_k;
end

geom = struct();
geom.q = q;
geom.rvec_in = vertcat(r_chunks{:});
geom.rout = vertcat(t_chunks{:});
geom.rin_body = rin_body;
geom.rout_body = rout_body;
geom.source_ind = source_ind;
geom.target_ind = target_ind;
geom.pairs = pairs;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Aii = Aii;
basis.use_fmm = use_fmm;
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

res = laplaceSingleLayerField(geom.rvec_in,geom.rout,lambda_all,basis.use_fmm);

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

function test_solve_laplace_1B
fprintf('--- solve_laplace_1B self-test ---\n');

opt = getLaplace2Dparams();
R = opt.alpha;
q = [0; 2*R+0.08*R; 6*R+1.5i*R];
g_body = [1; -0.7; 0.25];

[Q_it,lam_it,it_it,~,res_it] = solve_laplace_1B(q,g_body,0,1e-10,0,true);
[geom,~,~,~] = prepareLaplace1B(q,false);

rhs = zeros(length(geom.rout),1);
for k = 1:length(q)
    rhs(geom.target_ind{k}) = g_body(k);
end

A = lapSLPmat(geom.rvec_in,geom.rout);
[Yd,Ud] = getPseudoFactors(A,1e-14,0);
lam_dense = Yd*(Ud'*rhs);

Q_dense = zeros(length(q),1);
for k = 1:length(q)
    Q_dense(k) = sum(lam_dense(geom.source_ind{k}));
end

n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(length(q)*n_bound,1);
for k = 1:length(q)
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+R*(cos(tb)+1i*sin(tb));
end
rcheck_ext = buildExteriorPoints(q,R,600);

u_it_b = laplaceSingleLayerField(geom.rvec_in,rcheck_b,lam_it,false);
u_dense_b = laplaceSingleLayerField(geom.rvec_in,rcheck_b,lam_dense,false);

u_it_ext = laplaceSingleLayerField(geom.rvec_in,rcheck_ext,lam_it,false);
u_dense_ext = laplaceSingleLayerField(geom.rvec_in,rcheck_ext,lam_dense,false);

rel_Q = norm(Q_it-Q_dense,inf)/max(1,norm(Q_dense,inf));
rel_b = norm(u_it_b-u_dense_b,inf)/max(1,norm(u_dense_b,inf));
rel_ext = norm(u_it_ext-u_dense_ext,inf)/max(1,norm(u_dense_ext,inf));

fprintf('Iterative: it=%d, maxres=%.3e\n',it_it,res_it);
fprintf('Dense-vs-iterative rel diff Q      : %.3e\n',rel_Q);
fprintf('Dense-vs-iterative rel diff bndry  : %.3e\n',rel_b);
fprintf('Dense-vs-iterative rel diff ext    : %.3e\n',rel_ext);

end
