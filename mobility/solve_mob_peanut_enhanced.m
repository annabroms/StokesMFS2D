function [UW,sol] = solve_mob_peanut_enhanced(q,F,T,opt)
%SOLVE_MOB_PEANUT_ENHANCED Solve a 2D Stokes mobility problem with
%peanut-compressed pair corrections. The enhanced discretisation for each
%pair consists of Stokeslets only. 
%
% Syntax:
%   [UW,sol] = solve_mob_peanut_enhanced(q,F,T,opt)
%
% Inputs:
%   q          - Complex particle centers (P x 1).
%   F          - Net forces per particle (P x 2), columns are x y.
%   T          - Net torques per particle (P x 1).
%   opt        - Options struct.
%     Required geometry/discretisation fields:
%       rad           single particle radius (constant for all particles)
%       N_c,N_f       coarse/fine proxy source counts per particle
%       a_c           coarse collocation upsampling (M_c = ceil(a_c*N_c))
%       Rp_c,Rp_f     coarse/fine proxy radii, defined via tol_c by default
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
%       visualise_sol show diagnostic plots in postprocessing
%       surface_error_mode
%                     'rel' plots relative boundary errors, 'abs' plots
%                     absolute boundary errors
%       use_fmm       use FMM-accelerated Stokeslet evaluation in field
%                     evaluation paths that support it
%       cmap          if true, use coarse-to-coarse pair map without
%                     recovering fine sources (other than for evaluating the 
%                     field at or near the boundaries)
%       column_weight if true, scale least-squares matrix columns before
%                     SVD in the one-body, pair, and peanut factor builds
%       left_weight   if true, scale least-squares matrix rows by local
%                     arclength weights before SVD in the one-body, pair,
%                     and peanut factor builds
%
%
% Outputs:
%   UW         - Rigid-body motion vector [Ux1; Uy1; W1; ...; UxP; UyP; WP].
%   sol        - Struct with fields:
%                lambda_c, it, gmres_tol, rel_res, abs_res, resvec,
%                body_rel_res_max.
%
% Notes:
%   - Call with no inputs to run the built-in test.
%   - All Stokeslet sources contribute to the rigid body velocities
%     computed in UW.
%
% Anna Broms, Mar 2026

if nargin==0, test_solve_mob; 
    return; end

if nargin < 4 || ~isstruct(opt)
    error('solve_mob_peanut_enhanced requires q, F, T, and an options struct opt.');
end

q = q(:);
T = T(:);
P = numel(q);
assert(size(F,1)==P,'F must have one row per particle.');
assert(size(F,2)==2,'F must have two columns [Fx, Fy].');
assert(numel(T)==P,'T must have one entry per particle.');

visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = getOptField(opt,'gmres_tol',1e-10);
debug = logical(getOptField(opt,'debug',false));
surface_error_mode = getOptField(opt,'surface_error_mode','rel');
maxit = getOptField(opt,'maxit',800);
rad = getOptField(opt,'rad',1);
get_bndry_field = logical(getOptField(opt,'get_bndry_field',true));
get_precomp_time = logical(getOptField(opt,'get_precomp_time',false));
N_peanut = getOptField(opt,'N_peanut',400);
column_weight = logical(getOptField(opt,'column_weight',false));
left_weight = logical(getOptField(opt,'left_weight',false));
svd_opts = struct('column_weight',column_weight,'left_weight',left_weight);
precomp_time = struct('total',nan,'one_body',nan,'pair_setup',nan, ...
    'pair_basis',nan,'two_body_or_peanut',nan);
pair_setup_time = 0;
pair_basis_time = 0;
if N_peanut <= 0
    error('solve_mob_peanut_enhanced:InvalidNPeanut', ...
        ['solve_mob_peanut_enhanced requires opt.N_peanut > 0 ', ...
         'to build the peanut-compressed pair basis.']);
end
opt.N_peanut = N_peanut;

%% SET PARAMS
if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

%Set coarse and fine grid. 
%Play with N_c, N_f, a (a_f). 
N_c = getOptField(opt,'N_c',80);
N_f = getOptField(opt,'N_f',150);
a_c = getOptField(opt,'a_c',1.2);
tol_c = getOptField(opt,'tol_c',1e-10);

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rp_c = getOptField(opt,'Rp_c',rad*max([1-sep_c,0.01])); 
Rp_f = getOptField(opt,'Rp_f',rad*max([1-sep_f,0.01]));

%% Discretize
%Outer basic grid
tout_c_all = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c_all(1:end-1)';
rbase_out_c = rad*(cos(tout_c)+1i*sin(tout_c));
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rp_c*cos(tin)+Rp_c*1i*sin(tin); 

tin_f = linspace(0,2*pi,N_f+1);
tin_f = tin_f(1:end-1)';
rbase_in_f =  Rp_f*cos(tin_f)+Rp_f*1i*sin(tin_f);

if get_precomp_time
    pair_timer = tic;
end
[~, ~, ~, rimage_vec, refine,pairs] = getEnhancedGrid(q, opt);
if get_precomp_time
    pair_setup_time = toc(pair_timer);
    precomp_time.pair_setup = pair_setup_time;
end

rvec_in_c = [];
rout = [];

for k = 1:P
    rvec_in_c = [rvec_in_c; q(k)+rbase_in_c];
    rout = [rout; rbase_out_c+q(k)];
end

%% Prepare completion flow, one- and two-body factorisations
%get completion flow: evaluation of lambda0
[urhs,~] = getRecompletionFlow(rvec_in_c,rout,q,F,T); 

%Get one-body pseduo inverse blocks -- enough to do this for single body.
rimage_in = []; 
if get_precomp_time
    one_body_timer = tic;
end
[U,Y,Lc] = getSelfPseudoMobilityStokes( ...
    1,q,rbase_in_c,rbase_out_c,rimage_in,[0,ceil(a_c*N_c)],svd_opts);
if get_precomp_time
    precomp_time.one_body = toc(one_body_timer);
end

%Get pair basis
plot_grid = 0; %debug option: visualise each close pair 
opt.project_force = true;
opt.project = true;
opt.pair_basis_debug = plot_grid;
opt.rad = ones(P,1);
if get_precomp_time
    pair_timer = tic;
end
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU,pair_cache] = ...
    getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,Lc{1},rbase_out_c,svd_opts);
                     

% TODO: update visualisation:
%Visualise 1-body and pair-basis
%warning('Deactivate opt.precomp');
%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,[],[],N_c, N_f,a_c,a_f,rad)

Lc_pair = getILpair(Lc{1}); % Dense pair projector kept for fallback/comparison.
if get_precomp_time
    pair_basis_time = toc(pair_timer);
    precomp_time.pair_basis = pair_basis_time;
    precomp_time.two_body_or_peanut = pair_setup_time + pair_basis_time;
    precomp_time.total = precomp_time.one_body + precomp_time.two_body_or_peanut;
end

%% Construct check boundaries
% Create new grid points, for which the accuracy of the solution is
% to be evaluated. 
if get_bndry_field
    rcheck_b = [];
    n_bound = 803;
    t = linspace(0,2*pi,n_bound)';
    for k = 1:P
        rcheck_b = [rcheck_b; q(k)+rad*(cos(t)+1i*sin(t))];
    end
else
    rcheck_b = [];
    n_bound = 0;
end

%% Repackage to prepare for solve
% Shared input structs for matvec and transform calls.
% Use dot-assignment to keep cell fields intact and avoid struct-array expansion.
geom_solve = struct();
geom_solve.rbase_in_c = rbase_in_c;
geom_solve.rbase_in_f = rbase_in_f;
geom_solve.rvec_in = rvec_in_c;
geom_solve.rimage_vec = rimage_vec;
opt_solve = opt;
opt_solve.get_bndry_field = false;
geom_solve.opt = opt_solve;
geom_solve.rvec_out = rout;
geom_solve.rcheck = rout;
geom_solve.q = q;
geom_solve.pairs = pairs;
geom_solve.refine = refine; 

basis_mob = struct();
basis_mob.U = U;
basis_mob.Y = Y;
basis_mob.Lc = Lc{1};
basis_mob.Upf = UB_all;
basis_mob.Ypf = YB_all;
basis_mob.DC_all = UC_all;
basis_mob.YC_all = YC_all;
basis_mob.Cmap = Cmap;
basis_mob.Cmap_FU = Cmap_FU; 
basis_mob.Lc_pair = Lc_pair;
basis_mob.pair_cache = pair_cache;

geom_check = geom_solve;
geom_check.opt = opt;
if opt.cmap
    geom_solve.opt.get_bndry_field = 0;
end
geom_check.rcheck = rcheck_b;

%% Solve system

% Debug mode: build the matrix to check it out
if debug
    x = zeros(2*length(rout),1);
    tic
    ncols = 2*length(rout);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        uu = matvec_mob_peanut_enhanced(x,geom_solve,basis_mob);
        CC(:,k) = uu;
    end
    toc
    figure();
    imagesc(log10(abs(CC)))
    colorbar
    title([solver_name ': log_{10} |matvec system matrix|'],'interpreter','none')
    cc = skeel(CC);
    fprintf('Estimated condition number of system matrix: %1.3e \n',cc);
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
        quiver(real(rout),imag(rout),V_min(1:end/2,k),V_min(end/2+1:end,k));
        axis equal
    end
   


end

disp(' == Solving... == ');
[tau,it,resvec,real_res] = helsing_gmres(@(x) matvec_mob_peanut_enhanced(x,geom_solve,basis_mob),...
    urhs,2*size(rout,1),maxit,gmres_tol,opt.gmres_verbose,rout);

if visualise_sol
    figure()
    semilogy(resvec);
    xlabel('iteration number','interpreter','latex');
    ylabel('Estimated relative residual');
    axis tight
    grid on
    title('GMRES convergence with peanut compression, mobility', 'interpreter','latex')
end


if visualise_sol
    %check residual
    restot = (matvec_mob_peanut_enhanced(tau,geom_solve,basis_mob)-urhs)./urhs;
    figure()
    semilogy(abs(restot))
    title([solver_name ': Rel res at colloc points'],'interpreter','none')
    axis tight
end

disp(' == Postprocessing == ');
%% COMPUTE Rigid body motion
%And evaluate residual in new points rcheck_b

% Recover coarse and fine sources from data on the boundary / collocation grid.
geom_post = geom_solve;
if get_bndry_field
    geom_post = geom_check;
end

[lam_c_x, lam_c_nonpx,lam_self_x, lam_f_x,lam_c_y, ...
    lam_c_nonpy,lam_self_y,lam_f_y,~,rimage_k] = ...
    transform_mob_peanut_stokes(tau,geom_post,basis_mob);
lambda_c = [lam_c_x; lam_c_y];

%%% Get rigid body motion. 

%First due to all coarse sources
Kc = getKmat2D(rbase_in_c,0);
UW= zeros(3*P,1); 
for k= 1:P
    UW((k-1)*3+1:3*k) = -Kc'*[lam_c_nonpx((k-1)*N_c+1:k*N_c); lam_c_nonpy((k-1)*N_c+1:k*N_c)];
end

if opt.cmap
    for pair_it = 1:size(pairs,1)
        i = pairs(pair_it,1);
        p2 = pairs(pair_it,2);

        coarse_i = (i-1)*N_c+1:i*N_c;
        coarse_p2 = (p2-1)*N_c+1:p2*N_c;

        %Here the projected sources are used, as (I-L) is not yet
        %applied for Cmap_FU.
        rhs_pair = [lam_self_x(coarse_i); lam_self_x(coarse_p2); ...
                    lam_self_y(coarse_i); lam_self_y(coarse_p2)];

        % Determine rigid body motion for the pair, using ansatz
        if basis_mob.pair_cache.enabled
            pair = getStokesPairInstance(basis_mob.pair_cache,pair_it);
            rhs_pair = rotatePairOrderedStokesData(rhs_pair,N_c,pair.meta.phase_c,conj(pair.meta.rot));
            pair_vel = -pair.group.Cmap_FU*rhs_pair;
            vel_i = pair.meta.rot*(pair_vel(1) + 1i*pair_vel(2));
            vel_p2 = pair.meta.rot*(pair_vel(4) + 1i*pair_vel(5));
            pair_vel = [real(vel_i); imag(vel_i); pair_vel(3); ...
                        real(vel_p2); imag(vel_p2); pair_vel(6)];
        else
            pair_vel = -Cmap_FU{i,p2}*rhs_pair;
        end
        UW((i-1)*3+1:3*i) = UW((i-1)*3+1:3*i)+ pair_vel(1:3); 
        UW((p2-1)*3+1:3*p2) = UW((p2-1)*3+1:3*p2)+ pair_vel(4:6);

    end

else
    % Then, due to all pair sources (fine-body + fine-image).
    has_neigh = sort(unique(pairs(:)));
    for i = 1:length(has_neigh)
        k = has_neigh(i);
        rsrc_k = [rbase_in_f+q(k); rimage_k{k}];
        fx_k = lam_f_x{k};
        fy_k = lam_f_y{k};
        rel_k = rsrc_k - q(k);
    
        rbm_k = [sum(fx_k); ...
                 sum(fy_k); ...
                 sum((-imag(rel_k)).*fx_k + real(rel_k).*fy_k)];
        UW((k-1)*3+1:3*k) = UW((k-1)*3+1:3*k)-rbm_k;
    end
end

if get_bndry_field
    B = getKmat2D(rcheck_b(1:n_bound)-q(1),0); %same for all particles
    u_lhs = zeros(2*P*n_bound,1);
    for k = 1:P  
        res = B*UW(3*(k-1)+1:3*k);
        u_lhs((k-1)*n_bound+1:k*n_bound) = res(1:end/2);
        u_lhs(P*n_bound+(k-1)*n_bound+1:P*n_bound+k*n_bound) = res(end/2+1:end); 
    end
    geom_check.opt.self_correct = 0; 
    u_rhs = matvec_mob_peanut_enhanced(tau,geom_check,basis_mob);
    S_0 = getRecompletionFlow(rvec_in_c,rcheck_b,q,F,T); 
    u_rhs = u_rhs-S_0;

    disp('Surface residual')
    diff_vec = u_rhs-u_lhs;
    max_abs = max(abs(S_0(1:end/2)+1i*S_0(end/2+1:end)));
    res = abs(diff_vec(1:end/2)+1i*diff_vec(end/2+1:end));
    abs_res = max(res); 
    if max_abs > 0
        rel_vec = res/max_abs;
    else
        rel_vec = res;
    end
    rel_res = max(rel_vec);
    fprintf('Relative boundary error: %1.3e \n', rel_res);
    fprintf('Absolute boundary error: %1.3e \n', abs_res);
else
    rel_res = nan;
    abs_res = nan;
    rel_vec = [];
    res = [];
    body_rel_res_max = [];
    fprintf('Boundary field evaluation skipped (opt.get_bndry_field=0)\n');
end

if get_bndry_field
    rel_grid = reshape(rel_vec,n_bound,P);
    body_rel_res_max = max(rel_grid,[],1).';
end

 
if visualise_sol && get_bndry_field
    
    % visualise boundary velocities with some offset from boundary so that particles
    % visually don't overlap
    rvis = [];
    aa = 0.9;
    for k = 1:P
        rvis = [rvis; aa*(rcheck_b(n_bound*(k-1)+1:k*n_bound)-q(k))+q(k)];
    end

    if strcmp(surface_error_mode,'rel')
        err_plot = rel_vec;
        err_label = 'relative';
    else
        err_plot = res;
        err_label = 'absolute';
    end
    err_plot_log = log10(max(err_plot, eps));

    figure()
    scatter3(real(rvis),imag(rvis),err_plot_log,30,err_plot_log,'filled');
    hold on
    plot(real(rcheck_b),imag(rcheck_b),'k.')
    c = colorbar;
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    c.TickLabelInterpreter = 'latex';
    ylabel(c, sprintf('%s boundary error (log10)', err_label), 'interpreter','none','FontSize',18)
    c.FontSize = 18; 
    err_label_cap = [upper(err_label(1)) err_label(2:end)];
    title(sprintf('%s: %s boundary error', solver_name, err_label_cap), 'interpreter','none')


    %% Visualise source strengths
    figure()
    semilogy(abs(lambda_c))
    title('Coarse source strengths mobility, peanut compression')
    axis tight
    
end

sol = struct();
sol.lambda_c = lambda_c;
sol.it = it;
sol.gmres_unknowns = 2*length(rout);
sol.gmres_tol = gmres_tol;
sol.rel_res = rel_res;
sol.abs_res = abs_res;
sol.resvec = resvec;
sol.body_rel_res_max = body_rel_res_max;
sol.precomp_time = precomp_time;
sol.pair_precomp_stats = pair_cache.stats;


end


function test_solve_mob

close all;

%% Set geometry and data
q = [0; 2.001; 2.001i]; %center coordinates

delta = 1e-3;
%delta = 0.5; 
P = 20;
q = 0:2+delta:(P-1)*(2+delta);
% P = 4; 
% q = [0; 2+delta; 7; 9+delta];
%P = 20; 
side = 2 + delta;               % neighbor center distance
R = side / (2*sin(pi/P));         % ring radius
%q = R * exp(1i * (0:P-1).' * (2*pi/P));

%q(1) = 8;
%q = q+5; 
%q = q-q(1);´

rng(8);
%q = grow_cluster(P,delta,2);
%q = [0; 2+delta]*1i;
rings = 1;
q = hexagonal_lattice(delta,rings,1);
%q = [q(4);q(end-1:end)];

P = length(q);
F = [real(q) imag(q)];
T = zeros(size(q));
F = randn(P,2);
T = randn(P,1);
rad = ones(size(q));
%F = [1 0; -1 0];

%F = [1 0; 0 0; 0 1; -1 0]; %forces on the particles
%T = [1; 1; 1; -1]; %torques on the particles
%rad = [1; 1; 1; 1]; 

F = F-mean(F); %zero total force


delta_pair = 0.2;
%delta_pair = 2+0.4;
visualise = 1; % one-body solve specific
%% Solve
%[UW1,lambdahat,it1,gmres_tol, err1] = solve_mob_2B_images(q,F,T,rad,delta_pair,visualise);

%compare to a solution with image enhancement
N_peanut = 400;
%[UW1, lambdahat1,it1,gmres_tol, rel1, abs1] = solve_mob_2B_images(q,F,T,rad,delta_pair,visualise);
gmres_tol = 1e-8;
images = 1; 
lr = 0; % long-range precond
debug = 0;
%[UW1,lambda_mob,it1,gmres_tol,err1] = solve_mob_1B(q,F,T,rad,images, lr, visualise);
N_c = 150;
N_f = 150;
N_c = 60;
N_f = 60;
opt = get2Dparams(P,N_c,N_f);
opt.delta_pair = delta_pair;
opt.N_peanut = N_peanut;
opt.visualise_sol = 1;
opt.visualise_grid = 1; 
opt.gmres_tol = gmres_tol;
opt.debug = debug;
opt.surface_error_mode = 'rel';
opt.rotation_mode = 'oversampled_fft';
opt.reuse_pair_basis_by_sep = 1;
opt.cmap = 1; % coarse to coarse compression?
opt.self_correct = 1; % create identiy matrix for a pair by utilising known rhs in pair problem
opt.use_dense = 1; % use stored matrices for evaluation of Stokeslet on single body / pair
opt.pair_basis_debug = 0; 
% opt.beta = 0.5;
% opt.Nclust = 200;
[UWp,solp] = solve_mob_peanut_enhanced(q,F,T,opt); 
% opt.reuse_pair_basis_by_sep = 0;
% [UWp2,solp2] = solve_mob_peanut_enhanced(q,F,T,opt); 
% [UW3,lambdahat3,it3,gmres_tol, rel3, abs3] = solve_mob_peanut_images(q,F,T,rad,delta_pair,N_peanut,visualise,opt.gmres_tol,1);
% [UW2B,sol_2B] = solve_mob_2B_enhanced(q,F,T,opt);
% alignfigs(3);

end
