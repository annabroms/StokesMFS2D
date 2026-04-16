function [v_body,sol] = solve_elast_2B(q,Q_body,opt)
%SOLVE_ELAST_2B Solve exterior Laplace elastance problem (known charges, unknown voltages) with 2-body preconditioning.
%
% Syntax:
%   [v_body,sol] = solve_elast_2B(q,Q_body,opt)
%
% Inputs:
%   q          - Complex particle centers (P x 1).
%   Q_body     - Prescribed net charge per body (P x 1).
%   opt        - Options struct (see getLaplace2Dparams.m)
%     Required fields:
%       rad           physical particle radius
%       N_c,N_f       coarse/fine proxy point counts
%       a_c,a_f       coarse/fine collocation upsampling factors
%       Rp_c,Rp_f     coarse/fine proxy radii
%       delta_pair    pair-detection threshold
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
%       use_fmm       use fmm2d (of flatiron) for Laplace field evals
%       RAM_check     estimate/report RAM usage for precomp, solve, and
%                     postprocessing using memorygraph
%       visualise_sol plot postprocessing diagnostics
%       body_plot_font_size font size for the bodywise scalar plots
%
% Outputs:
%   v_body    - Recovered constant boundary values per body (P x 1).
%   sol       - Struct with fields:
%               lambda_all : stacked source strengths used in field evaluation
%               it         : GMRES iteration count
%               gmres_tol  : GMRES tolerance used
%               maxres     : max relative equipotential residual
%               resvec     : GMRES convergence history
%               ram_estimate : raw-byte RAM summary for the solver phases
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

if nargin < 3 || ~isstruct(opt)
    error('solve_elast_2B requires q, Q_body, and an options struct opt.');
end

[ram_check,ram_cleanup] = startRamCheck(opt,mfilename); %#ok<NASGU>

visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = getOptField(opt,'gmres_tol',1e-7);
debug = logical(getOptField(opt,'debug',false));
use_fmm = logical(getOptField(opt,'use_fmm',true));
body_plot_font_size = getOptField(opt,'body_plot_font_size',14);
get_precomp_time = logical(getOptField(opt,'get_precomp_time',false));
get_solve_time = logical(getOptField(opt,'get_solve_time',true));

q = q(:);
Q_body = Q_body(:);
P = numel(q);
assert(numel(Q_body)==P,'Q_body must have one entry per particle.');
precomp_time = struct('total',nan,'one_body',nan,'pair_setup',nan, ...
    'pair_basis',nan,'two_body_or_peanut',nan);

maxit = 800;

if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

rad = getOptField(opt,'rad',2);
N_c = getOptField(opt,'N_c',80);
N_f = getOptField(opt,'N_f',150);
a_c = getOptField(opt,'a_c',1.2);
a_f = getOptField(opt,'a_f',1.2);

tol_c = 1e-10;
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = getOptField(opt,'Rp_c',rad*max([1-sep_c,0.01]));
Rp_f = getOptField(opt,'Rp_f',rad*max([1-sep_f,0.01]));

opt.project_charge = true;

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

tout_f = linspace(0,2*pi,ceil(a_f*N_f)+1)';
tout_f = tout_f(1:end-1);
rout_base_f = rad*(cos(tout_f)+1i*sin(tout_f));

rout = zeros(P*nout,1);
rvec_in_c = zeros(P*N_c,1);
coarse_source_ind = cell(P,1);
for k = 1:P
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
    idx = (k-1)*N_c+1:k*N_c;
    coarse_source_ind{k} = idx;
    rvec_in_c(idx) = q(k)+rbase_in_c;
end

if get_precomp_time
    pair_setup_timer = tic;
end
[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
if get_precomp_time
    precomp_time.pair_setup = toc(pair_setup_timer);
end

%% Get 1- and 2-body basis 
if get_precomp_time
    pair_basis_timer = tic;
end
[Upf,Ypf,~,~,~,~,pair_cache] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f, ...
    rout_base_f,rbase_out_c,rimage_vec,refine,pairs,opt);
if get_precomp_time
    precomp_time.pair_basis = toc(pair_basis_timer);
    precomp_time.two_body_or_peanut = precomp_time.pair_setup + precomp_time.pair_basis;
end
if get_precomp_time
    one_body_timer = tic;
end
[UU,YY] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout],true);
if get_precomp_time
    precomp_time.one_body = toc(one_body_timer);
    precomp_time.total = precomp_time.one_body + precomp_time.two_body_or_peanut;
end

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
    title([solver_name ': log_{10}|matvec system matrix|'],'interpreter','none')
    [V,D] = eig(CC);
    D = diag(D);
    figure()
    plot(real(D),imag(D),'+')
    xlabel('Re \lambda')
    ylabel('Im \lambda')
    title([solver_name ': eigenvalues of matvec system matrix'],'interpreter','none')
end

ram_check = markRamCheckPhase(ram_check,'precomp_end');

disp(' == Solving... == ');
solve_time_token = manageSolveTimeMeasurement('start',get_solve_time);
solve_time_cleanup = onCleanup(@() manageSolveTimeMeasurement('reset'));
[tau,it,resvec,~] = helsing_gmres(@(x) matvec_lap_2B_enhanced(x,geom,basis,rout), ...
    u_rhs,length(rout),maxit,gmres_tol,opt.gmres_verbose,rout);
solve_time = manageSolveTimeMeasurement('finish',solve_time_token);
solve_time_cleanup = [];
ram_check = markRamCheckPhase(ram_check,'solve_end');

if visualise_sol
    figure(); semilogy(resvec)
    xlabel('iteration number','interpreter','latex');
    ylabel('Estimated relative residual');
    axis tight
    grid on
    title('GMRES convergence elastance 2B','interpreter','latex')
end

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

if visualise_sol
    figure();
    plot(u_b); hold on;
    plot(v_true);
    title('Boundary potential and per-body means (Laplace elastance 2B)')

    plotBodyScalars(q,rad,Q_body,v_body,body_plot_font_size, ...
        {'Given body charges','Computed body voltages'}, ...
        {'Charge','Voltage'}, ...
        'Bodywise charges and voltages');
end

sol = struct();
sol.lambda_all = lambda_all;
sol.it = it;
sol.gmres_tol = gmres_tol;
sol.maxres = maxres;
sol.resvec = resvec;
sol.precomp_time = precomp_time;
sol.pair_precomp_stats = pair_cache.stats;
sol.solve_time = solve_time;
sol.ram_estimate = finishRamCheck(ram_check);

end

function test_solve_elast_2B
fprintf('--- solve_elast_2B self-test ---\n');
close all;
run_two_way = false;


% Set geometry and data
rng(3);
R = 2;
delta = 1e-3;
P = 2;
q = grow_cluster(P,delta,2,R);
q = [0; (2+delta)*R]*1i*exp(1i*pi/4);
%hexagonal lattice
% x = R+R*delta/2;
% y = sqrt((2*R+delta*R)^2-(R+R*delta/2)^2);
% q = [0; 2*R+delta*R; x+1i*y];
mode = 4; 
q = hexagonal_lattice(delta,mode,R);
P = length(q); 
Q_body = randn(P,1);
Q_body = Q_body-mean(Q_body); %zero total charge

% Set params
opt = getLaplace2Dparams(P,R);
opt.rad = R;
opt.visualise_sol = 1;
opt.visualise_grid = 1; 
opt.debug = 0;
opt.gmres_tol = 1e-7;
opt.use_fmm = true;
opt.gmres_verbose = 0;
opt.reuse_pair_basis_by_sep=true;
opt.shared_sep_tol = 1e-2*delta*max(1,opt.rad);


% Solve
[v2,sol2] = solve_elast_2B(q,Q_body,opt);
opt_1B = opt;
opt_1B.visualise_sol = 0;
[v1,sol1] = solve_elast_1B(q,Q_body,opt_1B);
it2 = sol2.it;
res2 = sol2.maxres;
it1 = sol1.it;
res1 = sol1.maxres;

fprintf('2B: it=%d, maxres=%.3e\n',it2,res2);
fprintf('1B: it=%d, maxres=%.3e\n',it1,res1);
fprintf('Rel diff in v_body (2B vs 1B): %.3e\n',norm(v2-v1)/max(1,norm(v1)));

if run_two_way
    v_ref = randn(P,1);
    opt_tw = opt;
    opt_tw.visualise_sol = 0;
    opt_tw.debug = 0;
    opt_tw.gmres_tol = 1e-10;
    [Q_cap,~] = solve_cap_2B(q,v_ref,opt_tw);
    [v_back,~] = solve_elast_2B(q,Q_cap,opt_tw);
    rel_two = norm(v_back-v_ref,inf)/max(1,norm(v_ref,inf));
    fprintf('Two-way rel diff in v_body      : %.3e\n',rel_two);
end

alignfigs; 

end
