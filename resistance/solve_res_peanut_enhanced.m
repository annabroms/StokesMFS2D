function [FT,sol] = solve_res_peanut_enhanced(q,U,W,opt)
%SOLVE_PRECOND_PEANUT_ENHANCED Solves a 2D Stokes resistance problem with circular
%particles using a 2-body preconditioned MFS formulation. A
%fine grid enhanced with shielding Stokeslets near image points is used locally for every
%close pair of particles to resolve challenging interactions. This fine grid is
%compressed so that only coarse sources and collocation points are needed to 
%solve the problem iteratively, effectively preconditioning the system.
%
% Syntax:
%   [FT,sol] = solve_res_peanut_enhanced(q,U,W,opt)
%
% Inputs:
%   q          - Vector of length P, complex-valued center coordinates for the particles
%   U          - Px2 matrix of translational velocities (columns: x and y components)
%   W          - Px1 column vector of angular velocities
%   opt        - Options struct.
%     Required geometry/discretisation fields:
%       rad           single particle radius (constant for all particles)
%       N_c,N_f       coarse/fine proxy source counts per particle
%       a_c           coarse collocation upsampling (M_c = ceil(a_c*N_c))
%       Rp_c,Rp_f     coarse/fine proxy radii defined via tol_c by default
%       delta_pair    close-pair threshold used for enhanced discretisation
%       Nclust        total Chebyshev nodes on each enclosing ellipse used
%                     to extract shielding arcs for close pairs
%       beta          ellipse shape parameter for the shielding-node construction

%     Solver-control fields:
%       gmres_tol     GMRES tolerance
%       gmres_verbose GMRES print level:
%                     0 = silent
%                     1 = final convergence summary only
%                     2 = per-iteration residual history + final summary
%       maxit         maximum GMRES iterations
%       debug         build/plot/investigate system matrix corresponding to
%                     matvec.
%       use_fmm       use FMM-accelerated Stokeslet evaluation in field
%                     evaluation paths that support it
%       get_bndry_field
%                     if true, evaluate and report boundary velocity residuals
%       cmap          if true, use coarse-to-coarse pair map for FT updates
%       use_big_sparse
%                     if true, use solve-grid sparse close-pair correction
%                     maps in the GMRES matvec
%       column_weight if true, scale least-squares matrix columns before
%                     SVD in the one-body, pair, and peanut factor builds
%       left_weight   if true, scale least-squares matrix rows by local
%                     arclength weights before SVD in the one-body, pair,
%                     and peanut factor builds
%       RAM_check     estimate/report RAM usage for precomp, solve, and
%                     postprocessing using memorygraph
%       visualise_sol show diagnostic plots in postprocessing
%
% Outputs:
%   FT         - 3P×1 vector of computed net forces and torques 
%   sol        - Struct with fields:
%                lambda_proxy, it, gmres_tol, rel_res, resvec, ram_estimate.
%
% Description:
%   The FMM is used for Stokeslet evaluation. No other source types are
%   needed thanks to compression of the local fine grid.
%
% Notes:
%   - Aims to test an MFS generalisation of the idea presented by Cheng-Greengard (1998)
%
% See also:
%   solve_res_1B              - 1-body preconditioned resistance solver
%   solve_res_2B_enhanced     - 2-body preconditioner without peanut compression
%   solve_mob_peanut_enhanced - 2-body preconditioned mobility solver
%                               with peanut compression
%   matvec_res_peanut_enhanced - matrix-vector action used by GMRES
%
% To test: Call without arguments.
%
% Anna Broms, Mar 2026

if nargin==0, test_solve_res; 
    return; end

if nargin < 4 || ~isstruct(opt)
    error('solve_res_peanut_enhanced requires q, U, W, and an options struct opt.');
end

[ram_check,~] = startRamCheck(opt,mfilename); 

q = q(:);
W = W(:);
P = numel(q);
assert(size(U,1)==P,'U must have one row per particle.');
assert(size(U,2)==2,'U must have two columns [Ux, Uy].');
assert(numel(W)==P,'W must have one entry per particle.');

opt.resistance = 1; 
visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = opt.gmres_tol;
debug = opt.debug; 
maxit = opt.maxit;
rad = getOptField(opt,'rad',1);
get_precomp_time = opt.get_precomp_time;
get_solve_time = opt.get_solve_time;
column_weight = opt.column_weight;
left_weight = opt.left_weight;
svd_opts = struct('column_weight',column_weight,'left_weight',left_weight);
use_big_sparse = opt.use_big_sparse;
precomp_time = struct('total',nan,'one_body',nan,'pair_setup',nan, ...
    'pair_basis',nan,'big_sparse',nan,'two_body_or_peanut',nan);

if use_big_sparse
    if ~logical(getOptField(opt,'cmap',false))
        error('solve_res_peanut_enhanced:BigSparseUnsupported', ...
            'opt.use_big_sparse=1 requires opt.cmap=1.');
    end
    if ~logical(getOptField(opt,'self_correct',false))
        error('solve_res_peanut_enhanced:BigSparseUnsupported', ...
            'opt.use_big_sparse=1 requires opt.self_correct=1.');
    end
    if getOptField(opt,'N_peanut',0) <= 0
        error('solve_res_peanut_enhanced:BigSparseUnsupported', ...
            'opt.use_big_sparse=1 requires opt.N_peanut > 0.');
    end
    if logical(getOptField(opt,'reuse_pair_basis_by_sep',false))
        fprintf(['opt.use_big_sparse=1: disabling opt.reuse_pair_basis_by_sep ', ...
            'for the resistance solve-grid sparse matvec.\n']);
        opt.reuse_pair_basis_by_sep = false;
    end
end

%% SET PARAMS
if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

% Grid params
%Set coarse and fine grid. 
 
N_c = getOptField(opt,'N_c',150);
N_f = getOptField(opt,'N_f',150);
a_c = getOptField(opt,'a_c',1.2);
tol_c = getOptField(opt,'tol_c',1e-9);

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rp_c = getOptField(opt,'Rp_c',rad*max([1-sep_c,0.01]));
Rp_f = getOptField(opt,'Rp_f',rad*max([1-sep_f,0.01]));

%% CREATE GRID
%Outer basic grid
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';
rbase_out_c = cos(tout_c)+1i*sin(tout_c);

tin_f = linspace(0,2*pi,N_f+1);
tin_f = tin_f(1:end-1)'; 
rbase_in_f = Rp_f*cos(tin_f)+Rp_f*1i*sin(tin_f);
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rp_c*cos(tin)+Rp_c*1i*sin(tin); 
rvec_in_c = [];
for k = 1:P
    rvec_in_c = [rvec_in_c; q(k)+rbase_in_c];
end

% get discretization
rout = []; 
for k = 1:P
    rout = [rout; rbase_out_c+q(k)];
end
 
if get_precomp_time
    pair_setup_timer = tic;
end
[~, ~, ~, rimage_vec, refine,pairs] = getEnhancedGrid(q, opt);
if get_precomp_time
    precomp_time.pair_setup = toc(pair_setup_timer);
end


%Get pair basis 
opt.project_force = false;
if get_precomp_time
    pair_basis_timer = tic;
end
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU,pair_cache] = ...
    getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,[],rbase_out_c,svd_opts);
if get_precomp_time
    precomp_time.pair_basis = toc(pair_basis_timer);
    precomp_time.two_body_or_peanut = precomp_time.pair_setup + precomp_time.pair_basis;
end

%Get one-body pseduo inverse blocks -- enough to do this for single body.
if get_precomp_time
    one_body_timer = tic;
end
[UU,YY] = getSelfPseudo(1,rbase_in_c,rbase_out_c,[],[],[0,numel(rbase_out_c)],0,svd_opts);
if get_precomp_time
    precomp_time.one_body = toc(one_body_timer);
    precomp_time.total = precomp_time.one_body + precomp_time.two_body_or_peanut;
end

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.rimage_vec = rimage_vec;
geom.opt = opt;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.rvec_in = rvec_in_c;

basis = struct();
basis.U = UU;
basis.Y = YY;
basis.Upf = UB_all;
basis.Ypf = YB_all;
basis.DC_all = UC_all;
basis.YC_all = YC_all;
basis.Cmap = Cmap;
basis.Cmap_FU = Cmap_FU;
basis.pair_cache = pair_cache;
basis.Nii = stokSLPmat(rbase_in_c,rbase_out_c,1); 

%% Build solve-grid sparse pair maps when requested.
big_sparse_stats = initResBigSparseSolveStats( ...
    use_big_sparse,size(pairs,1));
if use_big_sparse
    [big_sparse,big_sparse_stats] = buildResPeanutBigSparseStokes( ...
        geom,basis);
    basis.big_sparse = big_sparse;
    geom.opt.use_big_sparse = true;
    if get_precomp_time
        precomp_time.big_sparse = big_sparse_stats.build_time;
        precomp_time.two_body_or_peanut = precomp_time.pair_setup + ...
            precomp_time.pair_basis + precomp_time.big_sparse;
        precomp_time.total = precomp_time.one_body + ...
            precomp_time.two_body_or_peanut;
    end
end

%% Construct rhs

%Set boundary condition for MFS. Evaluate known velocities on the coarse
%grid of collocation points. 
foutx = [];
fouty = [];
for k = 1:P
   rhs_x = @(x) U(k,1)-W(k)*(imag(x)-imag(q(k)));
   rhs_y = @(x) U(k,2)+W(k)*(real(x-q(k)));
   foutx = [foutx; rhs_x(rout((k-1)*ceil(a_c*N_c)+1:k*ceil(a_c*N_c)))];
   fouty = [fouty; rhs_y(rout((k-1)*ceil(a_c*N_c)+1:k*ceil(a_c*N_c)))];  
end

fout = [foutx; fouty];

%% Select the GMRES matvec.
if use_big_sparse
    matvec_handle = @(x) matvec_peanut_big_sparse(x,geom,basis);
else
    matvec_handle = @(x) matvec_res_peanut_enhanced(x,geom,basis);
end

%% Solve system

% Build the matrix to inspect conditioning/eigenvalues if requested.
if debug
    x = zeros(2*length(rout),1);
    tic
    ncols = 2*length(rout);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        uu = matvec_handle(x);
        CC(:,k) = uu;
    end
    toc
    figure(14);
    clf; 
    imagesc(log10(abs(CC)))
    colorbar
    title([solver_name ': log_{10} |matvec system matrix|'],'interpreter','none')
    cc = skeel(CC);
    fprintf('Estimated condition number of system matrix: %1.3e \n',cc);
    figure();
    [V,D] = eig(CC);
    D = diag(D); 
    plot(real(D),imag(D),'b+')
    title([solver_name ': eigenvalues of matvec system matrix'],'interpreter','none')


    num_min = 4;
    [e,ind] = mink(real(D),num_min);
    V_min = V(:,ind);
    figure()   
    for k = 1:num_min
        subplot(num_min,1,k);
        quiver(real(rout),imag(rout),V_min(1:end/2,k),V_min(end/2+1:end,k));
        axis equal
    end
end

ram_check = markRamCheckPhase(ram_check,'precomp_end');

disp(' == Solving... == ');
solve_time_token = manageSolveTimeMeasurement('start',get_solve_time);
solve_time_cleanup = onCleanup(@() manageSolveTimeMeasurement('reset'));
[tau,it,resvec,real_res] = helsing_gmres(matvec_handle, ...
    fout,2*size(rout,1),maxit,gmres_tol,opt.gmres_verbose,rout);
solve_time = manageSolveTimeMeasurement('finish',solve_time_token);
solve_time_cleanup = [];
ram_check = markRamCheckPhase(ram_check,'solve_end');

if visualise_sol
    figure()
    semilogy(resvec);
    xlabel('iteration number','interpreter','latex');
    ylabel('Estimated relative residual');
    axis tight
    grid on
    title('GMRES convergence with peanut compression, resistance', 'interpreter','latex')
end

if visualise_sol
    %check residual
    restot = matvec_handle(tau)-fout;
    figure()
    semilogy(abs(restot))
    title('Res at colloc points, peanut resistance')
    axis tight
end

disp(' == Postprocessing == ');
%% Reconstruct sources and optionally evaluate boundary velocities.

if opt.get_bndry_field
    % Create new grid points, for which the accuracy of the solution is
    % to be evaluated. 
    rcheck_b = [];
    n_bound = 803;
    t = linspace(0,2*pi,n_bound)';
    for k = 1:P
        rcheck_b = [rcheck_b; q(k)+cos(t)+1i*sin(t)];
    end

    % Boundary postprocessing intentionally uses the local contact loop.
    % No sparse boundary matrix is built for these off-solve-grid targets.
    geom_eval = geom;
    geom_eval.rcheck = rcheck_b;
    [lam_c_x, lam_self_x, lam_f_x,lam_c_y,lam_self_y,lam_f_y,u_corr,rimage_k] = ...
        transform_peanut_stokes(tau,geom_eval,basis);

    %% Do the evaluation of the flow in check points 
    ftest_b = getVelocityField(rvec_in_c,rcheck_b,lam_c_x,lam_c_y, ...
        logical(getOptField(opt,'use_fmm',true)));
    ftest_b = ftest_b+u_corr; % Apply corrections on pairs
    
    %% Compute the postprocessing boundary velocity residual.
    fbound_x = ftest_b(1:length(rcheck_b));
    fbound_y = ftest_b(length(rcheck_b)+1:end);
    
    fb_x = [];
    fb_y = [];
    
    for k = 1:P
        rhs_f = @(x) [U(k,1)-W(k)*(imag(x)-imag(q(k))); U(k,2)+W(k)*(real(x-q(k)))]; 
        fb_true = rhs_f(rcheck_b(n_bound*(k-1)+1:n_bound*k));
        fb_x = [fb_x; fb_true(1:n_bound)];
        fb_y = [fb_y; fb_true(n_bound+1:end)];   
    end
    
    rel_res = max(sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./max(sqrt(fb_x.^2+fb_y.^2));
    fprintf(['Max boundary velocity relative residual on ', ...
        'postprocessing nodes %.3e\n'], rel_res);
else
    geom_eval = geom;
    if use_big_sparse
        [lam_c_x, lam_c_y,lam_self,~] = ...
            transform_res_peanut_big_sparse_stokes(tau,geom_eval,basis);
        lam_self_x = reshape(lam_self(1:N_c,:),[],1);
        lam_self_y = reshape(lam_self(N_c+1:2*N_c,:),[],1);
        lam_f_x = [];
        lam_f_y = [];
    elseif opt.cmap
        [lam_c_x, lam_self_x, ~,lam_c_y,lam_self_y,~,~,rimage_k] = ...
            transform_peanut_stokes(tau,geom_eval,basis);
        lam_f_x = [];
        lam_f_y = [];
    else
        [lam_c_x, lam_self_x, lam_f_x,lam_c_y,lam_self_y,lam_f_y,~,rimage_k] = ...
            transform_peanut_stokes(tau,geom_eval,basis);
    end
    rel_res = nan; 
end

%return compressed coarse sources
lambda_proxy = [lam_c_x; lam_c_y];

%% Determine forces and torques
%from coarse sources
K = getKmat2D(rbase_in_c,0);
FT = zeros(3*P,1); 
for k= 1:P
    FT((k-1)*3+1:3*k) = K'*[lam_self_x((k-1)*N_c+1:k*N_c); lam_self_y((k-1)*N_c+1:k*N_c)];
end

% Then, due to all pair sources (fine-body + fine-image).
if opt.cmap && use_big_sparse && ...
        isfield(basis,'big_sparse') && isfield(basis.big_sparse,'M_ft_corr')
    FT_pair = basis.big_sparse.M_ft_corr*[lam_self_x; lam_self_y];
    FT = addBigSparseForceTorqueCorrections( ...
        FT,basis.big_sparse.ft_scatter_rows,FT_pair);
elseif opt.cmap
    % Use precomputed contact maps.
    for pair_it = 1:size(pairs,1)
        i = pairs(pair_it,1);
        p2 = pairs(pair_it,2);

        coarse_i = (i-1)*N_c+1:i*N_c;
        coarse_p2 = (p2-1)*N_c+1:p2*N_c;
        rhs_pair = [lam_self_x(coarse_i); lam_self_x(coarse_p2); ...
                    lam_self_y(coarse_i); lam_self_y(coarse_p2)];

        % Sign matches the fine-source postprocessing path via lam_f_{x,y}.
        if basis.pair_cache.enabled
            pair = getStokesPairInstance(basis.pair_cache,pair_it);
            rhs_pair = rotatePairOrderedStokesData(rhs_pair,N_c,pair.meta.phase_c,conj(pair.meta.rot));
            FT_pair = pair.group.Cmap_FU*rhs_pair;
            fi = pair.meta.rot*(FT_pair(1) + 1i*FT_pair(2));
            fj = pair.meta.rot*(FT_pair(4) + 1i*FT_pair(5));
            FT_pair = [real(fi); imag(fi); FT_pair(3); ...
                       real(fj); imag(fj); FT_pair(6)];
        else
            FT_pair = Cmap_FU{i,p2}*rhs_pair; % [Fx_i; Fy_i; T_i; Fx_p2; Fy_p2; T_p2]
        end
        FT((i-1)*3+1:3*i) = FT((i-1)*3+1:3*i) + FT_pair(1:3);
        FT((p2-1)*3+1:3*p2) = FT((p2-1)*3+1:3*p2) + FT_pair(4:6);
    end
else
    has_neigh = sort(unique(pairs(:)));
    for i = 1:length(has_neigh)
        k = has_neigh(i);
        src_k = [rbase_in_f+q(k); rimage_k{k}];
        Kpair = getKmat2D(src_k,q(k));
        FT((k-1)*3+1:3*k) = FT((k-1)*3+1:3*k)+Kpair'*[lam_f_x{k}; lam_f_y{k}];
    end
end


if visualise_sol

    if opt.get_bndry_field
        figure()
        plot(ftest_b)
        hold on
        plot([fb_x; fb_y])
        title('Velocity rhs and lhs, peanut resistance')
    
        % Visualise the error for x and y together
        figure()
        m = max(sqrt(fb_x.^2+fb_y.^2));
        scatter3(real(rcheck_b),imag(rcheck_b),log10((sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./m),...
            30,log10((sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./m),'filled')
        c = colorbar;
        axis equal
        view(0,90)
        grid off
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        ylabel(c,'Relative error in $u$ on boundary (log10)','Interpreter','latex')
        set(c,'TickLabelInterpreter','latex')
        set(gca,'TickLabelInterpreter','latex')
        set(gcf,'color','w');
        axis off
        title('Resistance peanut')
    end

    %% Visualise density
    figure()
    subplot(1,2,1)
    semilogy(abs(lambda_proxy))
    axis tight
    hold on
    subplot(1,2,2)
    plot(lambda_proxy)
    axis tight
    sgtitle('Source strengths $\lambda$ with peanut compression, resistance', 'interpreter','latex')

   
end

sol = struct();
sol.lambda_proxy = lambda_proxy;
sol.it = it;
sol.gmres_unknowns = 2*size(rout,1);
sol.gmres_tol = gmres_tol;
sol.real_res = real_res;
sol.rel_res = rel_res;
sol.resvec = resvec;
sol.precomp_time = precomp_time;
sol.pair_precomp_stats = pair_cache.stats;
sol.solve_time = solve_time;
sol.big_sparse_stats = big_sparse_stats;
sol.ram_estimate = finishRamCheck(ram_check);

end

function stats = initResBigSparseSolveStats(requested,n_pairs)
stats = struct();
stats.requested = logical(requested);
stats.active = false;
stats.backend = 'global_block_sparse';
stats.reason = 'not_requested';
stats.n_pairs = n_pairs;
stats.N_c = 0;
stats.N_check = 0;
stats.used_pair_cache = false;
stats.source_correction = '';
stats.velocity_correction = '';
stats.u_corr_mode = '';
stats.combined_u_corr = false;
stats.force_torque_correction = '';
stats.local_pair_entries = 0;
stats.local_u_entries = 0;
stats.local_u_cross_entries = 0;
stats.local_u_peanut_entries = 0;
stats.local_ft_entries = 0;
stats.nnz_pair = 0;
stats.nnz_u = 0;
stats.nnz_u_cross = 0;
stats.nnz_u_peanut = 0;
stats.nnz_ft = 0;
stats.nnz_source_scatter = 0;
stats.nnz_velocity_scatter = 0;
stats.nnz_ft_scatter = 0;
stats.big_sparse_matrix_bytes = 0;
stats.big_sparse_auxiliary_bytes = 0;
stats.big_sparse_build_bytes = 0;
stats.big_sparse_peak_bytes = 0;
stats.build_time = 0;
if requested
    stats.reason = 'not_prepared';
end
end

function FT = addBigSparseForceTorqueCorrections(FT,rows,FT_pair)
pair_rows = 6;
n_pairs = numel(FT_pair)/pair_rows;
for row = 1:n_pairs
    idx = (row-1)*pair_rows+1:row*pair_rows;
    FT(rows(idx)) = FT(rows(idx)) + FT_pair(idx);
end
end


function test_solve_res

close all; 
test = 2; 
delta_pair = 0.2;
%delta_pair = 2+0.4; 
N_peanut = 400; 
N_c = 150; 
N_f = 150; 
visualise = 1; 
debug = 1;

opt = get2Dparams(1,N_c,N_f);
opt.delta_pair = delta_pair;
opt.N_peanut = N_peanut;
opt.visualise_sol = visualise;
opt.debug = debug;
opt.cmap = 1; 
opt.reuse_pair_basis_by_sep = 1;
opt.rotation_mode = 'oversampled_fft';
opt.rotation_oversample = 8;
opt.visualise_sol = visualise;
opt.gmres_tol = 1e-8;
opt.debug = debug;
%opt.Nclust = 200; 
%opt.beta = 0.5;

%% Check that parameters make sense
report = test_pair_corrections_stokes(opt,@solve_mob_peanut_enhanced,@solve_res_peanut_enhanced);

%% Setup test and solve
if test == 1
    delta = 0.001;
    q = [0; 2+delta; (2+delta)*1i]; %center coordinates
    
    %or, instead, three circles in triangle - this DOES NOT converge to a good
    %sol! 
    x = 1+delta/2;
    y = sqrt((2+delta)^2-(1+delta/2)^2);
    q = [0; 2+delta; x+1i*y];
    
    
    U = [1 0; 0 0; 0 1]; %translational velocities for particles
    W = [1; 1; 1]; %angular velocities
    % U = U*1e-5;
    % W = W*1e-5;
    rad = [1; 1; 1]; 
    visualise = 0; 

    q = [0; 2+delta];
    %q = q+5;
    U = U(1:2,:);
    W = W(1:2); 
    rad = rad(1:2); 


    opt.P = 3;

    
    [FT1,sol1] = solve_res_peanut_enhanced(q,U,W,opt);
    debug = 0; 
    [FT2,lambda2,it2,gmres_res2, err2] = solve_res_peanut_images(q,U,W,rad,delta_pair,N_peanut,visualise,0,gmres_tol,debug);
    opt.debug = debug;
    [FT3,sol3] = solve_res_2B_enhanced(q,U,W,opt);
    
    F = [FT2(1:3:end) FT2(2:3:end)];
    T = FT2(3:3:end); 
    [UW,lambdahat,it1,gmres_mob, err_mob] = solve_mob_peanut_images(q,F,T,rad,delta_pair,N_peanut,visualise);
    Ures = [U W]';
    
    
    % Compare to solution with pair corrections but without peanut compression
    visualise = 0; 
    [FTp,lambda,it2,gmres_tol, err2] = solve_res_2B_images(q,U,W,rad,delta_pair,0,visualise,gmres_tol,debug);
    
    str = sprintf('Two way error is %1.3e',norm(Ures(:)-UW));
    disp(str); 
    str = sprintf('Relative residual with 2-body precond: %1.2e vs with peanut compression: %1.2e\n Converging in %u resp %u iterations', ...
        err2,sol1.rel_res,it2,sol1.it);
    disp(str)
    
    alignfigs(4);
else

    rng(9);
    P = 5;
    delta = 0.001; %P = 5
    x = 1+delta/2;
    y = sqrt((2+delta)^2-(1+delta/2)^2);
    q = [0; 2+delta; x+1i*y];
   % P = length(q); 

    q = grow_cluster(P,delta,2);
    q = (2+delta)*(0:P-1);
    rings = 1;
    q = hexagonal_lattice(delta,rings,1);
    P = length(q);
   % q = [q; -6+1.5i; -2-4i]; P = P+2;
    %q = q([1,2,4],:); P = 3; 
    U = rand(P,2); W = rand(P,1); 
    %W = zeros(P,1); 
    
    opt.P = P; 


   % opt.rotation_mode = 'fft';
    [FT1p,sol1p] = solve_res_peanut_enhanced(q,U,W,opt);
    opt.reuse_pair_basis_by_sep = 0;
    [FT2p,sol2p] = solve_res_peanut_enhanced(q,U,W,opt);
    [FT2,sol2] = solve_res_2B_enhanced(q,U,W,opt);
    
    str = sprintf('Relative residual with peanut compression: %1.2e vs 2B preconditioner without compression: %1.2e\n Converging in %u resp %u iterations', ...
        sol1p.rel_res,sol2.rel_res,sol1p.it,sol2.it);
    disp(str)
    
    alignfigs;




end

end
