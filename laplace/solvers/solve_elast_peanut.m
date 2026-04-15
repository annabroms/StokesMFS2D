function [v_body,sol] = solve_elast_peanut(q,Q_body,opt)
%SOLVE_ELAST_PEANUT Solve exterior Laplace elastance problem (known charges, unknown voltages) with peanut-compressed 2B preconditioning.
%
% Syntax:
%   [v_body,sol] = solve_elast_peanut(q,Q_body,opt)
%
% Inputs:
%   q          - Complex particle centers (P x 1).
%   Q_body     - Prescribed net charge per body (P x 1).
%   opt        - Options struct (see getLaplace2Dparams.m).
%     Required fields:
%       rad           physical particle radius
%       N_c,N_f       coarse/fine proxy point counts
%       a_c,a_f       coarse/fine collocation upsampling factors
%       Rp_c,Rp_f     coarse/fine proxy radii
%       delta_pair    pair-detection threshold
%       N_peanut      peanut boundary node count
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
%       cmap          use compressed coarse-to-coarse map
%       get_bndry_field
%                     if true, reconstruct boundary fields/residuals in
%                     postprocessing
%       RAM_check     estimate/report RAM usage for precomp, solve, and
%                     postprocessing using memorygraph
%       visualise_sol plot postprocessing diagnostics
%       body_plot_font_size font size for the bodywise scalar plots
%
% Outputs:
%   v_body      - Recovered constant voltage values per body (P x 1).
%   sol         - Struct with fields:
%                 lambda_proxy : coarse proxy strengths after correction merge
%                 it           : GMRES iteration count
%                 gmres_tol    : GMRES tolerance used
%                 maxres       : max relative equipotential residual
%                                (NaN if opt.get_bndry_field = 0)
%                 resvec       : GMRES convergence history
%                 ram_estimate : raw-byte RAM summary for the solver phases
%
% Notes:
%   The radius parameter is chosen with rad ~= 1 to avoid unit logarithmic
%   capacity in 2D.
%
% To test: call without inputs.
%
% See also: solve_elast_1B, solve_elast_2B, solve_cap_peanut, ...
%   transform_lap_peanut, matvec_lap_peanut_enhanced.
%
% Anna Broms, Mar 2026

if nargin==0
    test_solve_elast_peanut;
    return
end

if nargin < 3 || ~isstruct(opt)
    error('solve_elast_peanut requires q, Q_body, and an options struct opt.');
end

[ram_check,ram_cleanup] = startRamCheck(opt,mfilename); %#ok<NASGU>

visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = getOptField(opt,'gmres_tol',1e-7);
debug = logical(getOptField(opt,'debug',false));
use_fmm = logical(getOptField(opt,'use_fmm',true));
gmres_verbose = getOptField(opt,'gmres_verbose',0);
get_bndry_field = logical(getOptField(opt,'get_bndry_field',true));
body_plot_font_size = getOptField(opt,'body_plot_font_size',14);
get_precomp_time = logical(getOptField(opt,'get_precomp_time',false));

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
opt.gmres_verbose = gmres_verbose;

N_c = getOptField(opt,'N_c',80); %coarse proxy sources per body
N_f = getOptField(opt,'N_f',150); %fine proxy sources per body (used to construct pair corrections only)
a_c = getOptField(opt,'a_c',1.2); %a_c = M_c/N_c, where M_c is the number of coarse collocation nodes per body
a_f = getOptField(opt,'a_f',1.2);

% Set radii for proxy points (see Stein & Barnett 2022 for discussion on how to choose these)
tol_c = 1e-10;
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = getOptField(opt,'Rp_c',rad*max([1-sep_c,0.01]));
Rp_f = getOptField(opt,'Rp_f',rad*max([1-sep_f,0.01]));

opt.project_charge = true; % always true for elastance, false for capacitance
opt_solve = opt;
opt_solve.get_bndry_field = false;

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

rvec_in_c = zeros(P*N_c,1);
rout = zeros(P*nout,1);
coarse_source_ind = cell(P,1);
for k = 1:P
    idx = (k-1)*N_c+1:k*N_c;
    coarse_source_ind{k} = idx;
    rvec_in_c(idx) = q(k)+rbase_in_c;
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
end

if get_precomp_time
    pair_setup_timer = tic;
end
[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
if get_precomp_time
    precomp_time.pair_setup = toc(pair_setup_timer);
end

%% Get 1- and 2-body basis functions
if get_precomp_time
    pair_basis_timer = tic;
end
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_QV,pair_cache] = ...
    getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f,rbase_out_c, ...
    rimage_vec,refine,pairs,opt);
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
geom.opt = opt_solve;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.rimage_vec = rimage_vec;
geom.rvec_in = rvec_in_c;
geom.pair_cache = pair_cache;
if opt.cmap
    geom.opt.get_bndry_field = 0;
end

basis = struct();
basis.U = UU;
basis.Y = YY;
basis.Upf = UB_all;
basis.Ypf = YB_all;
basis.DC_all = UC_all;
basis.YC_all = YC_all;
basis.Cmap = Cmap;
basis.Cmap_QV = Cmap_QV;
basis.pair_cache = pair_cache;
basis.Nii = lapSLPmat(rbase_in_c,rbase_out_c);



%% Get rhs based on "completion flow"
[lambda0_c,u_rhs] = getChargeCompletionFlowLaplace(rvec_in_c,rout,coarse_source_ind,Q_body,use_fmm);

%% Solve
% Build system matrix via matvec with columns of the identity matrix as input, for debugging
if debug
    x = zeros(length(rout),1);
    CC = zeros(length(rout));
    ncols = length(rout);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvec_lap_peanut_enhanced(x,geom,basis);
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
[tau,it,resvec,~] = helsing_gmres(@(x) matvec_lap_peanut_enhanced(x,geom,basis), ...
    u_rhs,length(rout),maxit,gmres_tol,opt.gmres_verbose,rout);
ram_check = markRamCheckPhase(ram_check,'solve_end');

if visualise_sol
    figure(); semilogy(resvec)
    xlabel('iteration number','interpreter','latex');
    ylabel('Estimated relative residual');
    axis tight
    grid on
    title('GMRES convergence elastance peanut','interpreter','latex')
end

disp(' == Postprocessing == ');
%% Postprocess 
if get_bndry_field
    n_bound = 803;
    tb = linspace(0,2*pi,n_bound+1)';
    tb = tb(1:end-1);
    rcheck_b = zeros(P*n_bound,1);
    for k = 1:P
        rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+rad*(cos(tb)+1i*sin(tb));
    end

    geom_eval = geom;
    geom_eval.opt = opt;
    geom_eval.rcheck = rcheck_b;
    geom_eval.opt.get_bndry_field = 1; 
else
    geom_eval = geom;
end

[lam_c,~,~,~,u_corr,pair_qv_nonp,~,lam_self_nonp,lam_f_nonp,lam_e_nonp] = ...
    transform_lap_peanut(tau,geom_eval,basis);
lambda_proxy = lambda0_c+lam_c;

v_body = zeros(P,1);
for k = 1:P
    idx = (k-1)*N_c+1:k*N_c;
    if opt.cmap
        v_body(k) = -sum(lam_self_nonp(idx)) - pair_qv_nonp(k);
    else
        lambda_tot_k = [-lam_self_nonp(idx); -lam_f_nonp{k}; -lam_e_nonp{k}];
        v_body(k) = sum(lambda_tot_k);
    end
end

if get_bndry_field
    u_b1 = matvec_lap_peanut_enhanced(tau,geom_eval,basis);
    u_b2 = lapSLPfield(rvec_in_c,rcheck_b,lambda0_c,use_fmm); %+lam_c,use_fmm);
    %u_b = u_b+u_corr;
    u_b = u_b1+u_b2;

    v_true = zeros(P*n_bound,1);
    for k = 1:P
        v_true((k-1)*n_bound+1:k*n_bound) = v_body(k);
    end
    maxres = max(abs(u_b-v_true))/max(1,max(abs(v_true)));
    fprintf('Max relative equipotential residual at new nodes %.3e\n',maxres);
else
    u_b = [];
    v_true = [];
    maxres = nan;
    fprintf('Boundary field evaluation skipped (opt.get_bndry_field=0)\n');
end

if visualise_sol
    if get_bndry_field
        figure();
        plot(u_b); hold on;
        plot(v_true);
        title('Boundary potential and per-body means (Laplace elastance peanut)')
    end

    figure();
    semilogy(abs(lambda_proxy))
    title('Coarse proxy strengths (Laplace elastance peanut)')

    plotBodyScalars(q,rad,Q_body,v_body,body_plot_font_size, ...
        {'Given body charges','Computed body voltages'}, ...
        {'Charge','Voltage'}, ...
        'Bodywise charges and voltages');
end

sol = struct();
sol.lambda_proxy = lambda_proxy;
sol.it = it;
sol.gmres_tol = gmres_tol;
sol.maxres = maxres;
sol.resvec = resvec;
sol.precomp_time = precomp_time;
sol.pair_precomp_stats = pair_cache.stats;
sol.ram_estimate = finishRamCheck(ram_check);

end

function test_solve_elast_peanut
fprintf('--- solve_elast_peanut self-test ---\n');

close all; 
run_two_way = true;

rng(9);
%% Set geometry and data
R = 2; 
delta = 1e-2;
P = 20; 
q = grow_cluster(P,delta,2,R);
%q = [0; R*(2+delta)]*exp(1i*pi/13);
Q_body = randn(P,1);
Q_body = Q_body-mean(Q_body); %zero total charge

%% Tune parameters
N_c = 61;
N_f = 150;
opt = getLaplace2Dparams(P,R,N_c,N_f);
opt.visualise_sol = 0;
opt.debug = 0; % determine system matrix by using matvec with columns of the identity matrix as input
opt.use_fmm = true;
opt.gmres_verbose = 0; %no output from gmres
opt.delta_pair = 0.5; % largest distance where pair corrections are triggered
opt.N_peanut = 400; % nodes on peanut separation surface
opt.gmres_tol = 1e-8;
opt.compress_cmap = 0; %use low rank approximation of coarse-coarse map
opt.cmap_tol = 1e-10; 
opt.reuse_pair_basis_by_sep = 1;
opt.check_rotations = 0; 
opt.rotation_mode = 'fft'; 
opt.rotation_oversample = 8; 
opt.cmap = 1;
%% Check that parameters make sense
report = test_pair_corrections_laplace(opt,@solve_cap_peanut,@solve_elast_peanut);

%% Solve
[vp,solp] = solve_elast_peanut(q,Q_body,opt);
opt.reuse_pair_basis_by_sep = 0;
[vp2,solp2] = solve_elast_peanut(q,Q_body,opt);
opt_2B = opt;
opt_2B.visualise_sol = 0;
[v2,sol2] = solve_elast_2B(q,Q_body,opt_2B);
itp = solp.it;
resp = solp.maxres;
it2 = sol2.it;
res2 = sol2.maxres;

fprintf('Peanut: it=%d, maxres=%.3e\n',itp,resp);
fprintf('2B    : it=%d, maxres=%.3e\n',it2,res2);
fprintf('Rel diff in v_body (peanut vs 2B): %.3e\n',norm(vp-v2)/max(1,norm(v2)));

if run_two_way
    if opt.visualise_sol
        disp('Press any key to continue')
        pause();
    end
    v_ref = randn(P,1);
    opt_tw = opt;
    opt_tw.visualise_sol = 0;
    opt_tw.debug = 0;
    [Q_cap,~] = solve_cap_peanut(q,v_ref,opt_tw);
    [v_back,~] = solve_elast_peanut(q,Q_cap,opt_tw);
    rel_two = norm(v_back-v_ref,inf)/max(1,norm(v_ref,inf));
    fprintf('Two-way rel diff in v_body         : %.3e\n',rel_two);
end

end
