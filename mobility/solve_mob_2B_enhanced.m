function [UW,sol] = solve_mob_2B_enhanced(q,F,T,opt)
%SOLVE_MOB_PRECOND_ENHANCED Solves a 2D Stokes mobility problem with circular
%particles using a 2-body preconditioned recompleted MFS formulation. To resolve
% challenging close interactions, a fine 2-body BVP is solved for fine
% source strengths on each close pair of particles. The fine sources,
% correct the representation obtained from a coarse grid, effectively preconditioning the system.
%
% Syntax:
%   [UW,sol] = solve_mob_2B_enhanced(q,F,T,opt)
%
% Inputs:
%   q          - Vector of length P, complex-valued center coordinates for the particles
%   F          - Px2 matrix of net force vectors (columns: x and y components)
%   T          - Px1 column vector of torques acting on the particles
%   opt        - Options struct. Common fields:
%                delta_pair, visualise_sol, gmres_tol,opt.gmres_verbose,
%                surface_error_mode, use_fmm, N_c, N_f, a_c, a_f, tol_c,
%                RAM_check.
%       debug    build/plot/investigate system matrix corresponding to
%                matvec.
%
% Outputs:
%   UW         - 3P×1 vector of computed rigid-body motion (RBM) velocities
%   sol        - Struct with fields:
%                lambda, it, gmres_tol, rel_res, abs_res, resvec,
%                ram_estimate.
%
% Description:
%   This function applies a 2-body preconditioner (using pair corrections via local fine BVPs)
%   along with image systems to improve solution accuracy for near-contact particle configurations.
%   The FMM is used for Stokeslet evaluation, while other source types use direct summation.
%
% Notes:
%   - Intended to match solve_2D_precond_images, but applied to a mobility problem rather than resistance.
%
% See also:
%   solve_mob_1B              - 1-body preconditioned mobility solver
%   solve_mob_peanut_images - 2-body preconditioner with peanut
%                              compression (in this version of the algorithm, only coarse sources
%                              are needed at the solve stage).
%   solve_res_2B_enhanced   - 2-body preconditioned resistance solver
%
% To test: Call without arguments.
%
% Anna Broms, Mar 2026

if nargin==0, test_solve_mob; 
    return; end

if nargin < 4 || ~isstruct(opt)
    error('solve_mob_2B_enhanced requires q, F, T, and an options struct opt.');
end

[ram_check,ram_cleanup] = startRamCheck(opt,mfilename); %#ok<NASGU>

q = q(:);
T = T(:);
P = numel(q);
assert(size(F,1)==P,'F must have one row per particle.');
assert(size(F,2)==2,'F must have two columns [Fx, Fy].');
assert(numel(T)==P,'T must have one entry per particle.');

visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = getOptField(opt,'gmres_tol',1e-10);
debug = logical(getOptField(opt,'debug',false));
surface_error_mode = getOptField(opt,'surface_error_mode','abs');
gmres_verbose = getOptField(opt,'gmres_verbose',0);
maxit = getOptField(opt,'maxit',800);
get_precomp_time = logical(getOptField(opt,'get_precomp_time',false));
opt.rad = 1;
precomp_time = struct('total',nan,'one_body',nan,'pair_setup',nan, ...
    'pair_basis',nan,'two_body_or_peanut',nan);

surface_error_mode = lower(char(surface_error_mode));
if ~any(strcmp(surface_error_mode, {'abs','rel'}))
    error('surface_error_mode must be ''abs'' or ''rel''.')
end

%% SET PARAMS
if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

opt.project_force = true; % project out contribution to force / torque
opt.gmres_verbose = gmres_verbose;
%Params for coarse and fine grid. 
N_c = getOptField(opt,'N_c',150);
N_f = getOptField(opt,'N_f',150);
a_c = getOptField(opt,'a_c',1.2);
tol_c = getOptField(opt,'tol_c',1e-12);

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rp_c = getOptField(opt,'Rp_c',max([1-sep_c,0.01]));
Rp_f = getOptField(opt,'Rp_f',max([1-sep_f,0.01]));
opt.rad = 1;

%% CREATE GRID
%Outer basic grid
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';

rbase_out_c = cos(tout_c)+1i*sin(tout_c);
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rp_c*cos(tin)+Rp_c*1i*sin(tin);

tin_f = linspace(0,2*pi,N_f+1);
tin_f = tin_f(1:end-1)';
rbase_in_f =  Rp_f*cos(tin_f)+Rp_f*1i*sin(tin_f);
 
if get_precomp_time
    pair_setup_timer = tic;
end
[~, ~, ~, rimage_vec, refine,pairs] = getEnhancedGrid(q, opt);
if get_precomp_time
    precomp_time.pair_setup = toc(pair_setup_timer);
end

rin = [];
rout = [];

for k = 1:P
    rin = [rin; q(k)+rbase_in_c];
    rout = [rout; rbase_out_c+q(k)];
end

% % To illustrate for presentation / paper
% for k = 1:P
%     rin_f = [rin_f; q(k)+rbase_in_f]; %just for visuals
%     rout_f = [rout_f; q(k)+rbase_out_f]; %just for visuals
% end
% if visualise_sol
%     %Visualise the compression between coarse and fine grid
%     figure()
%     fill(real(rout_f(1:end/2)),imag(rout_f(1:end/2)),[0 0 0],'FaceAlpha',0.1,'EdgeColor','none')
%     hold on
%     fill(real(rout_f(end/2+1:end)),imag(rout_f(end/2+1:end)),[0 0 0],'FaceAlpha',0.1,'EdgeColor','none')
%     %set colors for presentation... 
%     red = [231,51,57]/256;
%     yellow = [251,186,0]/256;
%     blue = [25,105,189]/256;
%     newred = 0.7*red+0.3*yellow;
%     plot(real(rin),imag(rin),'k.');
%     hold on
%     %plot(real(rin_f),imag(rin_f),'.','Color',newred);
%     
%     
%     rf1 = refine{1,2};
%     rf2 = refine{2,1};
%     %plot(real(rout_f),imag(rout_f),'.','Color',blue)
%     %plot(cos(rf1),sin(rf1),'.','Color',blue)
%     %plot(cos(rf2)+real(q(2)),sin(rf2)+imag(q(2)),'.','Color',blue)
%     
%     plot(real(rout),imag(rout),'b.')
%     rimage1 = rimage_vec{2,1};
%     rimage2 = rimage_vec{1,2};
%     %plot(real(rimage1),imag(rimage1),'.','Color',newred);
%     %plot(real(rimage2),imag(rimage2),'.','Color',newred);
%     axis equal
% end

%% PREPARE PRECONDITIONING AND RHS

%get evaluation of lambda0, the completion sources, computed from known net forces and torques
[urhs,lambda_x,lambda_y] = getRecompletionFlow(rin,rout,q,F,T); 

% The pseudo inverse for a single particle is supposed to use the 
% coarse proxy grid + basic grid of collocation points only. We assume
% everyone has the same discretisation, so it's enough to do this once.
rimage_in = []; 
if get_precomp_time
    one_body_timer = tic;
end
[U,Y,Lc] = getSelfPseudoMobilityStokes(1,q,rbase_in_c,rbase_out_c,rimage_in,[0,ceil(a_c*N_c)]);
if get_precomp_time
    precomp_time.one_body = toc(one_body_timer);
end

%Get pair basis
opt.project = true;
opt.N_peanut = 0; %no peanut compression here!
if get_precomp_time
    pair_basis_timer = tic;
end
[Upf,Ypf,~,~,~,Cmap_FU,pair_cache] = ...
    getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,Lc{1},rbase_out_c);
if get_precomp_time
    precomp_time.pair_basis = toc(pair_basis_timer);
    precomp_time.two_body_or_peanut = precomp_time.pair_setup + precomp_time.pair_basis;
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

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Lc = Lc{1};
basis.Upf = Upf;
basis.Ypf = Ypf;
basis.Cmap_FU = Cmap_FU;
basis.pair_cache = pair_cache;


% Now, check pair basis up to the boundary. Is it nice and smooth?
%warning('Deactivate opt.precomp');
%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,Lc{1},Lf,N_c,N_f,a_c,a_f,rad)

%checkOneBasis(rbase_in_c,U,Y,Lc{1},N_c,a_c)
%% Construct check boundaries
% Create new grid points, in which the accuracy of the solution is
% to be evaluated. 
rcheck_b = [];
n_bound = 803;
t = linspace(0,2*pi,n_bound)';
for k = 1:P
    rcheck_b = [rcheck_b; q(k)+(cos(t)+1i*sin(t))];
end

%% SOLVE SYSTEM
% Build the matrix to check it out

if debug
    x = zeros(size(urhs));
    ncols = size(urhs,1);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        uu = matvec_mob_2B_enhanced(x,geom,basis);
        CC(:,k) = uu;
    end
    figure()
    imagesc(log10(abs(CC)))
    colorbar
    title([solver_name ': log_{10} |matvec system matrix|'],'interpreter','none')
    cc = skeel(CC);
    fprintf('Estimated condition number of system matrix: %1.3e \n',cc);
    get_nullspace = 0;
    get_eigs = 1; 

    if get_eigs

        [V,D] = eig(CC);
        D = diag(D); 
        figure()
        plot(real(D),imag(D),'+')
        xlabel('Re \lambda')
        ylabel('Im \lambda')
        title([solver_name ': eigenvalues of matvec system matrix'],'interpreter','none')
            
        [s,I] = mink(abs(D),3);

        if get_nullspace
            Vsmall = V(:,I).*s';
        
            [UU,S,VV] = svd(CC);
            SS = diag(S);
            plot_singular_values = false;
            if plot_singular_values
                figure()
                semilogy(SS)
                semilogy(SS,'+')
            end
            [s2,I2] = mink(SS,6);
        
            UUs = UU(:,I2);
            VVs = VV(:,I2); 
        
            Mc = round(a_c*N_c);
            t = linspace(0,2*pi,Mc+1);
            t = t(1:end-1)'; 
        
            %Visualise left and right singular vectors
            for i = 1:6
                for k = 1:P
                    Vpx = VVs((k-1)*Mc+1:k*Mc,i);
                    Vpy = VVs((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
                    Vp = abs(Vpx+1i*Vpy);
                    figure(33)
                    subplot(2,3,i)
                    scatter3(real(q(k))+cos(t),imag(q(k))+sin(t),Vp,40,Vp,'filled');
                    hold on
                    colorbar
                    view(0,90)
                    sgtitle('Right sing vec','interpreter','latex')
                    title(i)
                    axis off
                    axis square
        
                    figure(34)
                    subplot(2,3,i)
                    quiver(real(q(k))+cos(t),imag(q(k))+sin(t),Vpx,Vpy)
                    hold on
                    q2 = [q; q(1)];
                    plot(real(q2),imag(q2),'k-')
                    axis off
                    sgtitle('Right sing vec','interpreter','latex')
                    title(i)
                    axis square
        
        
                    Vpx = UUs((k-1)*Mc+1:k*Mc,i);
                    Vpy = UUs((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
                    Vp = abs(Vpx+1i*Vpy);
                    figure(35)
                    subplot(2,3,i)
                    scatter3(real(q(k))+cos(t),imag(q(k))+sin(t),Vp,40,Vp,'filled');
                    hold on
                    colorbar
                    view(0,90)
                    sgtitle('Left sing vec','interpreter','latex')
                    title(i)
                    axis off
                    axis square
        
                    figure(36)
                    subplot(2,3,i)
                    quiver(real(q(k))+cos(t),imag(q(k))+sin(t),Vpx,Vpy)
                    hold on
                    q2 = [q; q(1)];
                    plot(real(q2),imag(q2),'k-')
                    axis off
                    sgtitle('Left sing vec','interpreter','latex')
                    title(i)
                    axis square
                end
            end
                    
            %Visualise eigenvectors
        
            for i = 1:3
             
                for k = 1:P  
                    Vpx = Vsmall((k-1)*Mc+1:k*Mc,i);
                    Vpy = Vsmall((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
                    Vp = abs(Vpx+1i*Vpy);
                    figure(31)
                    subplot(1,3,i)
                    scatter3(real(q(k))+cos(t),imag(q(k))+sin(t),Vp,40,Vp,'filled');
                    hold on
                    colorbar
                    view(0,90)
                    sgtitle('Abs of eigvec','interpreter','latex')
                    axis off
                    axis square
        
        
                    figure(32)
                    subplot(1,3,i)
                    quiver(real(q(k))+cos(t),imag(q(k))+sin(t),Vpx,Vpy)
                    hold on
                    q2 = [q; q(1)];
                    plot(real(q2),imag(q2),'k-')
                    sgtitle('Eigvec','interpreter','latex')
                    axis off
                    axis square
                end
        
    
            end
        end

    end
    

end

ram_check = markRamCheckPhase(ram_check,'precomp_end');

disp(' == Solving... == ');
[tau,it,resvec,real_res] = helsing_gmres(@(x) matvec_mob_2B_enhanced(x,geom,basis),urhs,2*length(rout),maxit,gmres_tol,opt.gmres_verbose,rout);
ram_check = markRamCheckPhase(ram_check,'solve_end');
plot_gmres = visualise_sol;

%Modify to build with krylov preconditioning
%[tau, e2, precond] = precond_gmres(@(x) ...)
%it = length(e2); 

if plot_gmres
      figure()
%     semilogy(e2);
      semilogy(resvec);
      xlabel('iteration number','interpreter','latex');
      ylabel('Estimated relative residual');
      axis tight
      grid on
      title('GMRES convergence mobility, pair corrections', 'Interpreter','latex')
      u2 = matvec_mob_2B_enhanced(tau,geom,basis);
end 

disp(' == Postprocessing == ');
%% POSTPROCESS
% Recover the projected source representation and the nonprojected
% strengths needed for rigid-body postprocessing.

[rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y, ...
    tau_stokes_nonpx, tau_stokes_nonpy,tau_stokes_e_nonpx, tau_stokes_e_nonpy, rimage_k] = ...
    getMobPairTransformationStokes(tau,geom,basis);

lambda = [tau_stokes_x; tau_stokes_y]; %This are the fine source strengths

%% Get rigid body motion. 

%First due to all coarse sources
Kc = getKmat2D(rbase_in_c,0);
UW= zeros(3*P,1); 
for k= 1:P
    UW((k-1)*3+1:3*k) = -Kc'*[tau_stokes_nonpx((k-1)*N_c+1:k*N_c); tau_stokes_nonpy((k-1)*N_c+1:k*N_c)];
end

if logical(getOptField(opt,'cmap',false)) && (~isempty(basis.Cmap_FU) || basis.pair_cache.enabled)
    for row = 1:size(pairs,1)
        i = pairs(row,1);
        p2 = pairs(row,2);
        rhs_pair = [tau_stokes_x(coarse_ind{i}); tau_stokes_x(coarse_ind{p2}); ...
                    tau_stokes_y(coarse_ind{i}); tau_stokes_y(coarse_ind{p2})];
        if isfield(basis,'pair_cache') && basis.pair_cache.enabled
            pair = getStokesPairInstance(basis.pair_cache,row);
            rhs_pair = rotatePairOrderedStokesData(rhs_pair,opt.N_c,pair.meta.phase_c,conj(pair.meta.rot));
            pair_vel = -pair.group.Cmap_FU*rhs_pair;
            vel_i = pair.meta.rot*(pair_vel(1) + 1i*pair_vel(2));
            vel_p2 = pair.meta.rot*(pair_vel(4) + 1i*pair_vel(5));
            pair_vel = [real(vel_i); imag(vel_i); pair_vel(3); ...
                        real(vel_p2); imag(vel_p2); pair_vel(6)];
        else
            pair_vel = -basis.Cmap_FU{i,p2}*rhs_pair;
        end
        UW((i-1)*3+1:3*i) = UW((i-1)*3+1:3*i) + pair_vel(1:3);
        UW((p2-1)*3+1:3*p2) = UW((p2-1)*3+1:3*p2) + pair_vel(4:6);
    end
else
    has_neigh = sort(unique(pairs(:)));
    Kf = getKmat2D(rbase_in_f,0); 
    for i = 1:length(has_neigh)
        k = has_neigh(i); 
        UW((k-1)*3+1:3*k) = UW((k-1)*3+1:3*k)-Kf'*[tau_stokes_nonpx((k-1)*N_f+1+P*N_c:k*N_f+P*N_c); 
            tau_stokes_nonpy((k-1)*N_f+1+P*N_c:k*N_f+P*N_c)];
        Kim = getKmat2D(rimage_k{k},q(k));
        UW((k-1)*3+1:3*k) = UW((k-1)*3+1:3*k)-Kim'*[tau_stokes_e_nonpx{k}; tau_stokes_e_nonpy{k}];
    end
end

%% CHECK RESIDUAL AT SURFACE
%Compute velocity at surface
B = getKmat2D(rcheck_b(1:n_bound)-q(1),0); %same for all particles
u_lhs = zeros(2*P*n_bound,1);
for k = 1:P  
    res = B*UW(3*(k-1)+1:3*k);
    u_lhs((k-1)*n_bound+1:k*n_bound) = res(1:end/2);
    u_lhs(P*n_bound+(k-1)*n_bound+1:P*n_bound+k*n_bound) = res(end/2+1:end); 
end

% Using the recovered fine/coarse source representation directly.
u_rhs = getVelocityField(rvec_in,rcheck_b,tau_stokes_x,tau_stokes_y);

S_0 = getRecompletionFlow(rin,rcheck_b,q,F,T); 
u_rhs = u_rhs-S_0;  %Note! Sign here due to how we have defined the completion flow. 
                    %This is accordinng to the representation of the flow

%debug
% figure()
% subplot(1,2,1)
% plot(u_rhs);
% hold on
% plot(u_lhs);
% subplot(1,2,2)
% semilogy(abs(u_rhs-u_lhs));


disp('Surface residual')
diff_vec = u_rhs-u_lhs;
max_abs = max(abs(S_0(1:end/2)+1i*S_0(end/2+1:end)));
res = abs(diff_vec(1:end/2)+1i*diff_vec(end/2+1:end));
abs_res = max(res);
if max_abs > 0
    err_vec = res/max_abs;
else
    err_vec = res;
end
rel_res = max(err_vec);
fprintf('Relative boundary error: %.3e\n', rel_res);
fprintf('Absolute boundary error: %.3e\n', abs_res);


%Some visualisations
if visualise_sol

    % Visualise velocity error on different boundary to avoid that
    % particles visually overlap
    rvis = [];
    aa = 0.9;
    for k = 1:P
        rvis = [rvis; aa*(rcheck_b(n_bound*(k-1)+1:k*n_bound)-q(k))+q(k)];
    end

    if strcmp(surface_error_mode,'rel')
        err_plot = err_vec;
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
    ylabel(c, sprintf('%s boundary error (log10)', err_label), 'interpreter','none','FontSize',14)
    c.FontSize = 14;
    err_label_cap = [upper(err_label(1)) err_label(2:end)];
    title(sprintf('%s: %s boundary error', solver_name, err_label_cap), 'interpreter','none')


    %% Visualise source strengths
    figure()
    semilogy(abs(lambda))
    title('Source strengths mobility with pair corr')
    axis tight
    
end

sol = struct();
sol.lambda = lambda;
sol.it = it;
sol.gmres_tol = gmres_tol;
sol.rel_res = rel_res;
sol.abs_res = abs_res;
sol.resvec = resvec;
sol.precomp_time = precomp_time;
sol.pair_precomp_stats = pair_cache.stats;
sol.ram_estimate = finishRamCheck(ram_check);

end


function test_solve_mob

close all;
delta = 0.001;
rng(5);
%q = [0; 2+delta; (2+delta)*1i]; %center coordinates
q = [0; 2+delta]*exp(1i*pi/23);
%q = [0; 2+delta; 6; 8+delta;];
%or, instead, three cireles in triangle
% delta = 0.001; 
% x = 1+delta/2;
% y = sqrt((2+delta)^2-(1+delta/2)^2);
% q = [0; 2+delta; x+1i*y];

% F = [1 0; 0 0; 0 1]; %forces on the particles
% T = [1; 1; 1]; %torques on the particles
% rad = [1; 1; 1]; 
P = 5;
q = grow_cluster(P,delta,2); 

P = length(q); 
F = rand(P,2);
T = rand(P,1);
rad = ones(P,1);




%If only two particles
% q = 5+[0; 2+delta];
% F = F(1:2,:); 
% T = T(1:2); 
% rad = [1;1]; 

%F = F-mean(F); %zero total force

visualise = 0; 
images = 1; 
delta_pair = 0.2; 
lr= 0; 
%[UW1,lambda_1,it1,~,err1] = solve_mob_1B(q,F,T,rad,images, lr, visualise);

%compare to a solution with image enhancement
gmres_tol = 1e-8;
debug = 0; 
N_c = 60; 
opt = get2Dparams(length(q),N_c);
opt.delta_pair = delta_pair;
opt.visualise_sol = visualise;
opt.gmres_tol = gmres_tol;
opt.debug = debug;
opt.surface_error_mode = 'rel';
opt.visualise_grid = 0; 
opt.cmap = 1;
opt.reuse_pair_basis_by_sep = 1;
opt.visualise_grid = 1; 
opt.pair_basis_debug = 1;
opt.rotation_mode = 'oversampled_fft';
opt.rotation_mode = 'fft';
opt.rotation_oversample = 32; 
opt.project_coarse_in_reference_frame = 0; 

% is this a way to debug?
% opt.N_peanut = 0;%  
% opt.cmap = 0; 
% [UW2,sol2] = solve_mob_peanut_enhanced(q,F,T,opt);
close all;
[UW2,sol2] = solve_mob_2B_enhanced(q,F,T,opt);

str = sprintf('Relative residual with 1-body precond: %1.2e vs 2-body: %1.2e\n Converging in %u resp %u iterations', ...
    err1,sol2.rel_res,it1,sol2.it);
disp(str)
rel_diff_UW = norm(UW1(:)-UW2(:))/max(norm(UW1(:)),eps);
fprintf('Relative difference between UW1 and UW2: %1.2e\n', rel_diff_UW);

alignfigs;

end
