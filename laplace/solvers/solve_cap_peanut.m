function [Q,sol] = solve_cap_peanut(q,v_body,opt)
%SOLVE_CAP_PEANUT Solve exterior Dirichlet Laplace (capacitance:
%known voltages, unknown charges) with peanut-compressed 2B preconditioner.
%
% Syntax:
%   [Q,sol] = solve_cap_peanut(q,v_body,opt)
% Inputs:
%   q          - Complex particle centers (P x 1).
%   v_body     - Constant boundary values per body (P x 1).
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
%       cmap          use compressed coarse to coarse map
%       get_bndry_field
%                     if true, reconstruct boundary fields/residuals in
%                     postprocessing
%       RAM_check     estimate/report RAM usage for precomp, solve, and
%                     postprocessing using memorygraph
%       visualise_sol plot postprocessing diagnostics
%       body_plot_font_size
%                     base font size used in the voltage/charge subplot figure
%
% Outputs:
%   Q           - Per-body unweighted sums of source strengths.
%   sol         - Struct with fields:
%                 lambda_proxy : compressed coarse source strengths
%                 it           : GMRES iteration count
%                 gmres_tol    : GMRES tolerance used
%                 maxres       : max relative residual on independent boundary points
%                                (NaN if opt.get_bndry_field = 0)
%                 resvec       : GMRES convergence history
%                 ram_estimate : raw-byte RAM summary for the solver phases
%
% Notes:
%   The radius parameter should be chosen with rad ~= 1 to avoid unit logarithmic
%   capacity in 2D.
%
% To test: call without inputs.
%
% See also: solve_cap_1B, solve_cap_2B, solve_elast_peanut, ...
%   transform_lap_peanut, matvec_lap_peanut_enhanced.
%
% Anna Broms, Mar 2026

if nargin==0
    test_solve_cap_peanut;
    return
end

if nargin < 3 || ~isstruct(opt)
    error('solve_cap_peanut requires q, v_body, and an options struct opt.');
end

[ram_check,] = startRamCheck(opt,mfilename);

visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = getOptField(opt,'gmres_tol',1e-7);
debug = logical(getOptField(opt,'debug',false));
use_fmm = logical(getOptField(opt,'use_fmm',true));
gmres_verbose = getOptField(opt,'gmres_verbose',0);
get_bndry_field = logical(getOptField(opt,'get_bndry_field',true));
body_plot_font_size = getOptField(opt,'body_plot_font_size',14);
get_precomp_time = logical(getOptField(opt,'get_precomp_time',false));
get_solve_time = logical(getOptField(opt,'get_solve_time',true));
use_big_sparse_requested = logical(getOptField(opt,'use_big_sparse',false));
big_sparse_build_mode = '';
if use_big_sparse_requested
    if ~logical(getOptField(opt,'cmap',false))
        error('solve_cap_peanut:BigSparseUnsupported', ...
            'opt.use_big_sparse=1 requires opt.cmap=1.');
    end
    big_sparse_build_mode = lower(char(getOptField(opt, ...
        'lap_big_sparse_build_mode','auto')));
    switch big_sparse_build_mode
        case 'auto'
            if get_bndry_field
                big_sparse_build_mode = 'precomputed';
            else
                big_sparse_build_mode = 'streaming';
            end
        case {'precomputed','streaming'}
            % supported
        otherwise
            error('solve_cap_peanut:BadBigSparseBuildMode', ...
                ['opt.lap_big_sparse_build_mode must be ''auto'', ', ...
                 '''precomputed'', or ''streaming''.']);
    end
    if strcmp(big_sparse_build_mode,'streaming') && get_bndry_field
        error('solve_cap_peanut:StreamingBoundaryUnsupported', ...
            ['lap_big_sparse_build_mode=''streaming'' requires ', ...
             'opt.get_bndry_field=0. Use ''auto'' or ''precomputed'' ', ...
             'for boundary postprocessing.']);
    end
end
opt_solve = opt;
opt_solve.get_bndry_field = false;
if use_big_sparse_requested
    opt_solve.use_big_sparse = true;
    opt_solve.lap_big_sparse_build_mode = big_sparse_build_mode;
end

q = q(:);
v_body = v_body(:);
P = numel(q);
assert(numel(v_body)==P,'v_body must have one entry per particle.');
precomp_time = struct('total',nan,'one_body',nan,'pair_setup',nan, ...
    'pair_basis',nan,'big_sparse',nan,'two_body_or_peanut',nan);~


%revert single threaded settings
setenv('OMP_NUM_THREADS','');
setenv('MKL_NUM_THREADS','');
setenv('OPENBLAS_NUM_THREADS','');

maxNumCompThreads('automatic');

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

rvec_in_c = zeros(P*N_c,1);
rout = zeros(P*nout,1);
for k = 1:P
    rvec_in_c((k-1)*N_c+1:k*N_c) = q(k)+rbase_in_c;
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
end

if get_precomp_time
    pair_setup_timer = tic;
end
[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
if get_precomp_time
    precomp_time.pair_setup = toc(pair_setup_timer);
end

%% Basis factors
streaming_big_sparse = use_big_sparse_requested && ...
    strcmp(big_sparse_build_mode,'streaming');
if get_precomp_time
    pair_basis_timer = tic;
end
if streaming_big_sparse
    UB_all = [];
    YB_all = [];
    UC_all = [];
    YC_all = [];
    Cmap = [];
    Cmap_QV = [];
    pair_cache = initLaplacePairCache();
    pair_cache.enabled = logical(getOptField(opt,'reuse_pair_basis_by_sep',false));
    if pair_cache.enabled
        pair_cache.stats.branch = 'canonical_group';
    else
        pair_cache.stats.branch = 'per_pair';
    end
    pair_cache.stats.n_pairs = size(pairs,1);
else
    [UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_QV,pair_cache] = ...
        getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f,rbase_out_c, ...
        rimage_vec,refine,pairs,opt);
end
if get_precomp_time
    precomp_time.pair_basis = toc(pair_basis_timer);
    precomp_time.two_body_or_peanut = precomp_time.pair_setup + precomp_time.pair_basis;
end
if get_precomp_time
    one_body_timer = tic;
end
[UU,YY] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout]);
if get_precomp_time
    precomp_time.one_body = toc(one_body_timer);
    precomp_time.total = precomp_time.one_body + precomp_time.two_body_or_peanut;
end

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.rout_base_f = rout_base_f;
geom.refine = refine;
geom.opt = opt_solve;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.rimage_vec = rimage_vec;
geom.rvec_in = rvec_in_c;
geom.pair_cache = pair_cache;

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

big_sparse_stats = initLaplaceBigSparseSolveStats( ...
    use_big_sparse_requested,size(pairs,1));
if use_big_sparse_requested
    [big_sparse,big_sparse_stats,pair_cache] = ...
        buildLaplacePeanutBigSparse(geom,basis);
    basis.big_sparse = big_sparse;
    basis.pair_cache = pair_cache;
    geom.pair_cache = pair_cache;
    geom.opt.use_big_sparse = true;
    if get_precomp_time
        precomp_time.big_sparse = big_sparse_stats.build_time;
        if streaming_big_sparse
            precomp_time.pair_basis = big_sparse_stats.build_time;
            precomp_time.two_body_or_peanut = precomp_time.pair_setup + ...
                precomp_time.pair_basis;
        else
            precomp_time.two_body_or_peanut = precomp_time.pair_setup + ...
                precomp_time.pair_basis + precomp_time.big_sparse;
        end
        precomp_time.total = precomp_time.one_body + ...
            precomp_time.two_body_or_peanut;
    end
end

%% RHS
fout = zeros(P*nout,1);
for k = 1:P
    fout((k-1)*nout+1:k*nout) = v_body(k);
end

%% Solve
if use_big_sparse_requested
    matvec_handle = @(x) matvec_lap_peanut_big_sparse(x,geom,basis);
else
    matvec_handle = @(x) matvec_lap_peanut_enhanced(x,geom,basis);
end

if debug
    x = zeros(length(rout),1);
    CC = zeros(length(rout));
    ncols = length(rout);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvec_handle(x);
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

    num_min = 3;
    [e_min,ind] = mink(abs(D),num_min);
    V_min = V(:,ind);
    figure()   
    for k = 1:num_min
        subplot(num_min,1,k);
        scatter(real(rout),imag(rout),10,log10(abs(V_min(:,k))),"filled");
        colorbar
        axis equal
    end
end

single_threaded = opt.single_threaded; 
if single_threaded
    %not sure if all of this is needed...
    setenv('OMP_NUM_THREADS','1');      % OpenMP
    setenv('MKL_NUM_THREADS','1');      % MATLAB/Intel MKL
    setenv('OPENBLAS_NUM_THREADS','1'); % OpenBLAS

    maxNumCompThreads(1);               % MATLAB computational threads

end

ram_check = markRamCheckPhase(ram_check,'precomp_end');

disp(' == Solving... == ');
solve_time_token = manageSolveTimeMeasurement('start',get_solve_time);
solve_time_cleanup = onCleanup(@() manageSolveTimeMeasurement('reset'));
[tau,it,resvec,~] = helsing_gmres(matvec_handle, ...
    fout,length(rout),maxit,gmres_tol,opt.gmres_verbose,rout);
solve_time = manageSolveTimeMeasurement('finish',solve_time_token);
ram_check = markRamCheckPhase(ram_check,'solve_end');

if visualise_sol
    figure(); semilogy(resvec)
    xlabel('iteration number','interpreter','latex');
    ylabel('Estimated relative residual');
    axis tight
    grid on
    title('GMRES convergence capacitance peanut','interpreter','latex')
end

disp(' == Postprocessing == ');
%% Postprocess
if get_bndry_field
    n_bound = 803;
    tb = linspace(0,2*pi,n_bound+1)';
    tb = tb(1:end-1);
    rcheck_b = zeros(P*n_bound,1);
    for k = 1:P
        rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+R*(cos(tb)+1i*sin(tb));
    end

    geom_eval = geom;
    geom_eval.opt = opt;
    geom_eval.rcheck = rcheck_b;
else
    geom_eval = geom;
end

if use_big_sparse_requested && isequal(geom_eval.rcheck,geom_eval.rvec_out)
    [lam_c,~,u_corr,pair_qv_nonp,~,lam_self_nonp] = ...
        transform_lap_peanut_big_sparse(tau,geom_eval,basis);
else
    [lam_c,~,~,~,u_corr,pair_qv_nonp,~,lam_self_nonp, ...
        lam_f_nonp,lam_e_nonp] = transform_lap_peanut(tau,geom_eval,basis);
end

lambda_proxy = lam_c;

if get_bndry_field
    u_b = lapSLPfield(rvec_in_c,rcheck_b,lam_c,use_fmm) + u_corr;

    g_true = zeros(P*n_bound,1);
    for k = 1:P
        g_true((k-1)*n_bound+1:k*n_bound) = v_body(k);
    end

    maxres = max(abs(u_b-g_true))/max(1,max(abs(g_true)));
    fprintf('Max surface relative residual at new nodes %.3e\n',maxres);
else
    u_b = [];
    g_true = [];
    maxres = nan;
    fprintf('Boundary field evaluation skipped (opt.get_bndry_field=0)\n');
end

Q = zeros(P,1);
for k = 1:P
    idx = (k-1)*N_c+1:k*N_c;
    if opt.cmap
        Q(k) = sum(lam_self_nonp(idx)) + pair_qv_nonp(k);
    else
        Q(k) = sum(lam_self_nonp(idx)) + sum(lam_f_nonp{k}) + sum(lam_e_nonp{k});
    end
end

if visualise_sol
    if get_bndry_field
        figure();
        plot(u_b); hold on; plot(g_true)
        title('Boundary values: lhs vs rhs (capacitance peanut)', ...
            'Interpreter','latex')
        set(gca,'TickLabelInterpreter','latex')
    end

    figure();
    semilogy(abs(lambda_proxy))
    title('Compressed source strengths (capacitance peanut)', ...
        'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    axis tight

    plotBodyScalars(q,R,v_body,Q,body_plot_font_size);
end

sol = struct();
sol.lambda_proxy = lambda_proxy;
sol.it = it;
sol.gmres_tol = gmres_tol;
sol.maxres = maxres;
sol.resvec = resvec;
sol.precomp_time = precomp_time;
sol.pair_precomp_stats = pair_cache.stats;
sol.solve_time = solve_time;
sol.big_sparse_stats = big_sparse_stats;
sol.ram_estimate = finishRamCheck(ram_check);

end

function test_solve_cap_peanut
fprintf('--- solve_cap_peanut self-test ---\n');

close all; 
rng(8);

% Set geometry and data
R = 2;
P = 10;
delta = 1e-3;
q = grow_cluster(P,delta,2,R);
q = (2*R+delta*R)*(0:P-1);
offset = 1;
%q = q+(-offset/2*1i+offset*1i*rand(size(q)));

% Solve capacitance for hexagonal lattice
%rings = 8; 
%rings = 5; 
%q = hexagonal_lattice(delta,rings,R);
P = length(q); 
v_body = buildAlternatingVoltages(q,R);
v_body = rand(P,1);
check_multi_compress = 0; 

% Set parameters and settings
N_c = 60; 
N_f = 60; 
opt = getLaplace2Dparams(P,R,N_c,N_f);
%opt.delta_pair = 2*R+0.6;
%opt.delta_pair = 4*R+0.6; %not better
%opt.delta_pair = 0.4;
opt.Nclust = 100;
opt.N_peanut = 400;
opt.visualise_sol = 1;
opt.visualise_grid = 0; 
opt.gmres_tol = 1e-8;
opt.debug = 0;
opt.use_fmm = true;
opt.gmres_verbose = 0;
opt.compress_cmap = 0; %use low rank approximation of coarse-coarse map
opt.cmap_tol = 1e-8; % tolerance used in the low-rank compression
opt.reuse_pair_basis_by_sep = 0;
opt.get_bndry_field = 1; 

%tic
[Qp,solp] = solve_cap_peanut(q,v_body,opt);
%t_one = toc;

opt.visualise_sol = 0;
opt.reuse_pair_basis_by_sep = 0; 
if check_multi_compress
  %  tic 
    [Qp2,solp2] = solve_cap_peanut(q,v_body,opt);
   % t_recomp = toc;
    fprintf('Peanut solution times: one pair=%1.2d s, multipair=%1.2d s, iters: %u vs %u, soldiff: %1.2d\n',t_one,...
        t_recomp,solp.it,solp2.it,norm(norm(Qp-Qp2)/max(1,norm(Qp2))));
else
%    fprintf('Solved and postprocessed in %1.2d s\n',t_one);
end
[Q2,sol2] = solve_cap_2B(q,v_body,opt);
itp = solp.it;
resp = solp.maxres;
it2 = sol2.it;
res2 = sol2.maxres;

fprintf('Peanut: it=%d, maxres=%.3e\n',itp,resp);
fprintf('2B    : it=%d, maxres=%.3e\n',it2,res2);
fprintf('Rel diff in Q (peanut vs 2B): %.3e\n',norm(Qp-Q2)/max(1,norm(Q2)));

end
