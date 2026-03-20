function [Q,sol] = solve_cap_2B(q,v_body,opt)
%SOLVE_CAP_2B Solve exterior Dirichlet Laplace problem
%(capacitance: known voltages, unknown charges) with 2 body preconditioning
%(no compression).
%
% Syntax:
%   [Q,sol] = solve_cap_2B(q,v_body,opt)
% Inputs:
%   q          - Complex particle centers (P x 1).
%   v_body     - Constant boundary values per body (P x 1).
%   opt        - Options struct (see getLaplace2Dparams.m).
%     Required fields:
%       rad           physical particle radius
%       N_c,N_f       coarse/fine proxy point counts
%       a_c,a_f       coarse/fine collocation upsampling factors
%       Rp_c,Rp_f     coarse/fine proxy radii
%       delta_pair    pair compression threshold
%       Nclust        total Chebyshev nodes on each enclosing ellipse used
%                     to extract the shielding arc of enhancing sources for
%                     each close pair
%     Solver-control fields:
%       gmres_tol     GMRES tolerance
%       gmres_verbose GMRES print level:
%                     0 = silent, 1 = final summary only,
%                     2 = per-iteration estimated residuals + final summary
%       debug         build/plot/investigate system matrix corresponding to
%                     matvec.
%       visualise_sol plot postprocessing diagnostics
%       body_plot_font_size font size for the bodywise scalar plots
%       use_fmm       use fmm2d (of flatiron) for Laplace field evals
%
% Outputs:
%   Q          - Net charges per body, computed as per-body unweighted sums of source strengths.
%   sol        - Struct with fields:
%                lambda_all : stacked fine source strengths used in field evaluation
%                it         : GMRES iteration count
%                gmres_tol  : GMRES tolerance used
%                maxres     : max relative residual on independent boundary points
%                resvec     : GMRES convergence history
%
% Notes:
%   The radius parameter is chosen with rad ~= 1 to avoid unit logarithmic
%   capacity in 2D.
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

if nargin < 3 || ~isstruct(opt)
    error('solve_cap_2B requires q, v_body, and an options struct opt.');
end

visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = getOptField(opt,'gmres_tol',1e-7);
debug = logical(getOptField(opt,'debug',false));
use_fmm = logical(getOptField(opt,'use_fmm',true));
gmres_verbose = getOptField(opt,'gmres_verbose',0);
body_plot_font_size = getOptField(opt,'body_plot_font_size',14);

q = q(:);
v_body = v_body(:);
P = numel(q);

assert(numel(v_body)==P,'v_body must have one entry per particle.');

%% Parameters
maxit = 800;

if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

R = getOptField(opt,'rad',2);
opt.gmres_verbose = gmres_verbose;

N_c = getOptField(opt,'N_c',80);
N_f = getOptField(opt,'N_f',150);
a_c = getOptField(opt,'a_c',1.2);
a_f = getOptField(opt,'a_f',1.2);

tol_c = 1e-10;
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = getOptField(opt,'Rp_c',R*max([1-sep_c,0.01]));
Rp_f = getOptField(opt,'Rp_f',R*max([1-sep_f,0.01]));

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

tout_f = linspace(0,2*pi,ceil(a_f*N_f)+1)';
tout_f = tout_f(1:end-1);
rout_base_f = R*(cos(tout_f)+1i*sin(tout_f));

rout = zeros(P*nout,1);
for k = 1:P
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
end

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

%% Basis factors
[Upf,Ypf,~,~,~,~,pair_cache] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f, ...
    rout_base_f,rbase_out_c,rimage_vec,refine,pairs,opt);
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
geom.pair_cache = pair_cache;

basis = struct();
basis.U = UU;
basis.Y = YY;
basis.Upf = Upf;
basis.Ypf = Ypf;
basis.pair_cache = pair_cache;

%% RHS
fout = zeros(P*nout,1);
for k = 1:P
    fout((k-1)*nout+1:k*nout) = v_body(k);
end

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
    title([solver_name ': log_{10}|matvec system matrix|'],'interpreter','none')
    [V,D] = eig(CC);
    D = diag(D);
    figure()
    plot(real(D),imag(D),'+')
    xlabel('Re \lambda')
    ylabel('Im \lambda')
    title([solver_name ': eigenvalues of matvec system matrix'],'interpreter','none')
end

disp(' == Solving... == ');
[tau,it,resvec,~] = helsing_gmres(@(x) matvec_lap_2B_enhanced(x,geom,basis,rout), ...
    fout,length(rout),maxit,gmres_tol,opt.gmres_verbose,rout);

figure(); semilogy(resvec)
xlabel('iteration number','interpreter','latex');
ylabel('Estimated relative residual');
axis tight
grid on
title('GMRES convergence capacitance 2B','interpreter','latex')

disp(' == Postprocessing == ');
%% Postprocess
[rvec_in,coarse_ind,lambda_all,lam_c,lam_f,lam_e] = getPairTransformationLaplace(tau,geom,basis);

n_bound = 803;
tb = linspace(0,2*pi,n_bound+1)';
tb = tb(1:end-1);
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+R*(cos(tb)+1i*sin(tb));
end

u_b = lapSLPfield(rvec_in,rcheck_b,lambda_all,use_fmm);
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

if visualise_sol
    figure();
    plot(u_b); hold on; plot(g_true)
    title('Boundary values: lhs vs rhs (Laplace 2B)')

    plotBodyScalars(q,R,v_body,Q,body_plot_font_size);
end

sol = struct();
sol.lambda_all = lambda_all;
sol.it = it;
sol.gmres_tol = gmres_tol;
sol.maxres = maxres;
sol.resvec = resvec;

end

function test_solve_cap_2B
fprintf('--- solve_cap_2B self-test ---\n');

close all; 


%% Set geometry and data
rng(8);
delta = 1e-2; 
R = 2; 
P = 6;
q = grow_cluster(P,delta,2,R);

%hexagonal lattice
%x = R+R*delta/2;
%y = sqrt((2*R+delta*R)^2-(R+R*delta/2)^2);
%q = [0; 2*R+delta*R; x+1i*y];
% mode = 5; 
mode = 1;
q = hexagonal_lattice(delta,mode,R);
q =  [0; R*(2+delta);10+[0;R*(2+delta)]*1i];
% 
P = length(q); 
v_body = rand(P,1); 

%% Set params and settings
N_c = 40;
opt = getLaplace2Dparams(P,R,N_c);
opt.visualise_sol = 1;
opt.visualise_grid = 1; 
opt.gmres_tol = 1e-12;
opt.debug = 1;
opt.use_fmm = true;
opt.gmres_verbose = 0;
opt.Nclust = 100;

[Q2,sol2] = solve_cap_2B(q,v_body,opt);
[Q1,sol1] = solve_cap_1B(q,v_body,opt);
it2 = sol2.it;
res2 = sol2.maxres;
it1 = sol1.it;
res1 = sol1.maxres;

fprintf('2B: it=%d, maxres=%.3e\n',it2,res2);
fprintf('1B: it=%d, maxres=%.3e\n',it1,res1);
fprintf('Rel diff in Q (2B vs 1B): %.3e\n',norm(Q2-Q1)/max(1,norm(Q1)));

end
