function [UW,lambdahat,it,gmres_tol,rel_res,abs_res] = solve_mob_1B(q,F,T,rads,image,lr,visualise,gmres_tol,debug,surface_error_mode,gmres_verbose)
%SOLVE_MOB_1B Solves a 2D Stokes mobility problem with circular particles
% using 1-body preconditioned recompleted MFS, with a combination of source
% types at enhancing nodes
%
% Syntax:
%   [UW, lambdahat, it, gmres_tol, rel_res, abs_res] = solve_mob_1B(q, F, T, rads, image, visualise)
%
% Inputs:
%   q         - Vector of length P, complex-valued center coordinates for the particles
%   F         - Px2 matrix of net force vectors (columns: x and y components)
%   T         - Px1 column vector of torques acting on the particles
%   rads      - Px1 vector of particle radii
%   image     - Logical flag (true/false): use image system for higher accuracy
%   lr        - 0,1,2,..: sets use of long-range preconditioning. If false, no
%               lr precond. lr > 1 determines the coarse space.
%   visualise - Logical flag: plot the configuration and solution details
%   gmres_tol - Optional GMRES tolerance (default 1e-10)
%   debug     - Optional logical flag: build/draw dense matrix CC and its
%               eigenvalues (default false)
%   surface_error_mode - Optional boundary-error plot mode: 'abs' (default)
%               or 'rel'
%
% Outputs:
%   UW         - 3P×1 vector of computed rigid-body motion (RBM) velocities
%   lambdahat  - Solution vector of source strengths
%   it         - Number of GMRES iterations required using GMRES with low stagnation control (J. Helsing).
%   gmres_tol  - Set GMRES tolerance
%   rel_res    - Maximum relative residual in a test (non-collocation) set of boundary nodes
%   abs_res    - Maximum absolute residual in a test (non-collocation) set of boundary nodes
%
% See also:
%   solve_mob_2B_images - 2-body preconditioning enhanced with images for mobility
%   solve_mob_peanut_images - 2-body preconditioning with peanut compression for mobility
%   solve_res_1B             - 1-body preconditioned resistance formulation
%
% To test, call without arguments
%
% Anna Broms, April 2, 2025


if nargin==0, test_solve_mob; 
    return; end

P = length(q);
if nargin < 4 || isempty(rads), rads = ones(P,1); end
if nargin < 5 || isempty(image), image = 1; end
if nargin < 6 || isempty(lr), lr = 0; end
if nargin < 7 || isempty(visualise), visualise = 0; end
if nargin < 8 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 9 || isempty(debug), debug = false; end
if nargin < 10 || isempty(surface_error_mode), surface_error_mode = 'abs'; end
if nargin < 11 || isempty(gmres_verbose), gmres_verbose = 0; end
surface_error_mode = lower(char(surface_error_mode));
if ~any(strcmp(surface_error_mode, {'abs','rel'}))
    error('surface_error_mode must be ''abs'' or ''rel''.')
end


%% Checks

assert(size(T,1)==P,'Wrong size of torque vector')
assert(size(F,1)==P,'Wrong size of force vector')
assert(size(F,2)==2,'Wrong size of force vector, should contain x y coordinates')

if sum(rads)~=P
    warning('Double check radii') %in the design of image and collocation points, it's assumed all radii are 1.
end


%% PARAMS
%GMRES
maxit = 1600; 
solver_name = 'solve_mob_1B';

project_proxy = 1; %version of the algorithm project_proxy = 0 corresponds to the version where also image points are projected

opt = get2Dparams(); 
opt.gmres_verbose = gmres_verbose;

%Params to set Rp and grid resolution, based on heuristic relationship from
%QFS paper (Laplace).
a = opt.a_c; 
opt.image = image;
opt.lr = lr; 
opt.P = P; 

if image
    Nc = 60; %has been 60 all the time before
    %Nc = 200; 
    Nc = 100; 
    
   % Nc = 100; 
    tol = 1e-12;

    % Nc = 50; 
    %tol = 1e-8;
else
    Nc = 200;
    Nc = 1000;
    Nc = 100; 
    
    tol = 1e-14;
    %tol = 1e-12; 
    %tol = 1e-10; % tol here should probably be larger to allow for smaller coefficients, which in turn will improve on the GMRES convergence 
    %tol = 1e-14;
end

sep = (1/Nc)*log(1/tol);

Rp = max([1-sep,0.01]); %radius of proxy surface
%Rp = 0.8;
opt.Rp_c = Rp;
opt.N_c = Nc;
opt.pc = 0; %no pair correction


%solve with Stresslets + potential dipoles at image points
s = [0 0 1 1 0 0 0]; % s = [S R T D]
opt.s = s; 

%% GET GRIDS AND VISUALISE

%create grids
[rout,weights,rin,rimage,nimage,pair_points] = get2DImageGrid(q,rads,opt); 

 
if visualise
    figure()
    %Visualise particles
    plot(real(rin),imag(rin),'b.');
    hold on
    plot(real(rout),imag(rout),'k.');
    if size(rimage,1)
        plot(real(rimage),imag(rimage),'r.');
        quiver(real(rimage),imag(rimage),real(nimage),imag(nimage))
        legend('proxy sources','collocation points','images sources','Stresslet directions')
    else
        legend('proxy sources','collocation points')
    end
    title('Problem setup','interpreter','latex')
    axis equal
end

% Create new grid points for checking accuracy of residual
rcheck_b = [];
n_bound = 803;
for k = 1:length(q)
    t = linspace(0,2*pi,n_bound)';
    rcheck_b = [rcheck_b; q(k)+rads(k)*(cos(t)+1i*sin(t))];
end

%% PREPARE PRECONDITIONING AND RHS

%get evaluation of lambda0, the completion sources, computed from known net forces and torques
[u,lambda_x,lambda_y] = getRecompletionFlow(rin,rout,q,F,T); 

%Compute self-interaction blocks for each particle separately if images are in use.  
if ~isempty(rimage)
    [UU,Y,L] = getSelfPseudoMobility(P,q,rin,rout,rimage,nimage,pair_points,s,project_proxy);
else
    %Without images, all pseudoinverses are the same so it's enough to do
    %this for a single body (assuming every particle is discretised with
    %the same number of source and collocation points)
    [UU,Y,L,Lr] = getSelfPseudoMobility(1,0,rin(1:Nc)-q(1),rout(1:a*Nc)-q(1),rimage,nimage,pair_points,s,project_proxy);
end

%% Experiment with left preconditioner based on deflation
solve = 1;
if lr
    [Sinv,Zi,Yi,db] = get_long_range_precond_mob(q,rin,rout,L{1},Lr,opt);
    opt.db = db;
    tau_coarse = getCoarseSource(u,Sinv,Zi,Yi,db,P,Nc,a);

    disp('...Long range preconditioner constructed')

    % Check if this is enough for a solution to the system, by evaluating
    % and comparing lhs to rhs
    proj_tau = zeros(P*Nc*2,1); 
    rhs = zeros(P*Nc*a*2,1); 
    for k = 1:P
        tau_k = [tau_coarse(Nc*(k-1)+1:Nc*k); tau_coarse(P*Nc+Nc*(k-1)+1:P*Nc+Nc*k)];
        rhs_k = -Lr*tau_k;
        rhs(Nc*a*(k-1)+1:Nc*a*k) = rhs_k(1:end/2); 
        rhs(Nc*a*(k-1)+P*Nc*a+1:Nc*a*k+P*Nc*a) = rhs_k(end/2+1:end);
        proj_tau_k = tau_k-L{1}*tau_k;
        proj_tau(Nc*(k-1)+1:Nc*k) = proj_tau_k(1:end/2); 
        proj_tau(Nc*(k-1)+Nc*P+1:Nc*k+Nc*P) = proj_tau_k(end/2+1:end);
    end
    lhs = getVelocityField(rin,rout,proj_tau(1:end/2),proj_tau(end/2+1:end),[],[],[],[],[],[], []);
    if norm(lhs-rhs)<gmres_tol
        solve = 0; 
    end
end

% Build system matrix to inspect conditioning/eigenvalues if requested.
if debug && lr
    x = zeros(2*length(rout),1);
    tic
    for k = 1:2*length(rout)
        k
        x(:) = 0; 
        x(k) = 1; 
        uu = lr_matvec_2D_mobility(x,rin,rout,rimage,nimage,q,UU,Y,L,Lr,pair_points,s,Sinv,Nx,Ny,Mx,Zi,Yi,opt);
        %uu = matvec_mob_1B(x,rin,rout,rout,rimage,nimage,q,UU,Y,L,pair_points,s,1,project_proxy);
        CC(:,k) = uu;
    end
    toc
    figure(14);
    clf; 
    imagesc(log10(abs(CC)))
    colorbar
    title([solver_name ': log_{10} |CC|'],'interpreter','none')
    skeel(CC)

    %Check eigvals of system matrix
    figure()
    [~,D] = eig(CC);
    D = diag(D); 
    figure()
    plot(real(D),imag(D),'+')
    xlabel('Re \lambda')
    ylabel('Im \lambda')
    title([solver_name ': eigenvalues of CC'],'interpreter','none')
end


%% SOLVE
if opt.lr && solve
    disp('...Solving for fine component...')
    Pf = applyPmat_mob(u,rin,rout,L{1},Lr,Sinv,Zi,Yi,opt); 
    [mu_gmres,it,resvec,real_res] = helsing_gmres(@(x) lr_matvec_2D_mobility(x,rin,rout,rimage,nimage,q,UU,Y,L,Lr,pair_points,s,Sinv,Zi,Yi,opt),Pf,2*size(rout,1),maxit,gmres_tol,opt,rout);

elseif solve
    [mu_gmres,it,resvec,real_res] = helsing_gmres(@(x) matvec_mob_1B(x,rin,rout,rout,rimage,nimage,q,UU,Y,L,pair_points,s,1,project_proxy),u,2*size(rout,1),maxit,gmres_tol,opt,rout);
end
% Decay of residual with iteration number
figure()
semilogy(resvec)
title('GMRES convergence mobility, 1-body precond', 'Interpreter','latex')
xlabel('Iteration number','interpreter','latex')

%% GET LAMBDA
% We have now solved for unknows at the boundary. Want to go back to
% unknowns at the sources! 
N_small = size(rin,1)/P;
PM = length(rout);

[tau_stokes_x,tau_stokes_y,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y] = getTransformedDensity(mu_gmres,rimage,UU,Y,P,N_small,PM,pair_points,s,project_proxy);


if opt.lr
    tau_stokes = applyQmat_mob([tau_stokes_x; tau_stokes_y],rin,rout,L{1},Lr,Sinv,Zi,Yi,opt);
    tau_stokes_x = tau_stokes(1:end/2)+tau_coarse(1:end/2); 
    tau_stokes_y = tau_stokes(end/2+1:end)+tau_coarse(end/2+1:end); 
end

tau_proxy = [tau_stokes_x; tau_stokes_y];

tau_image = [tau_stress_x; tau_stress_y; tau_pot_x; tau_pot_y];
lambdahat = [tau_proxy; tau_image];

%% get rigid body motion. 
UW = zeros(3*P,1); 
if project_proxy
    K = getKmat2D(rin(1:N_small,:),q(1));
else
    start_ind = 0; 
end

for k = 1:P
    if project_proxy      
        UW(3*(k-1)+1:3*k) = -K'*[tau_stokes_x((k-1)*N_small+1:k*N_small); tau_stokes_y((k-1)*N_small+1:k*N_small)];
    else
        %If the projection is done also for the image sources w other
        %source types, it requires more bookkeeping to retrieve RBM
        %velocities... 
        rim = rimage(start_ind+1:start_ind+pair_points(k,1));

        K = getKmat2D([rin((k-1)*N_small+1:k*N_small); rim; rim],q(k));

        inds = start_ind+1:start_ind+pair_points(k,1);
        tau_x = tau_stokes_x((k-1)*N_small+1:k*N_small);
        tau_y = tau_stokes_y((k-1)*N_small+1:k*N_small);
        tau_st_x = tau_stress_x(inds);
        tau_st_y = tau_stress_y(inds);
        tau_p_x = tau_pot_x(inds);
        tau_p_y = tau_pot_y(inds);


        density = [tau_x; tau_st_x; tau_p_x; tau_y; tau_st_y; tau_p_y];
        UW(3*(k-1)+1:3*k) = -K'*density;
        start_ind = start_ind+pair_points(k,1);
    end

end

%% CHECK RESIDUAL AT SURFACE
%Compute velocity at surface
B = getKmat2D(rcheck_b(1:n_bound),q(1)); %same for all particles
for k = 1:P  
    res = B*UW(3*(k-1)+1:3*k);
    u_lhs((k-1)*n_bound+1:k*n_bound) = res(1:end/2); 
    u_lhs(P*n_bound+(k-1)*n_bound+1:P*n_bound+k*n_bound) = res(end/2+1:end); 
end

%Using representation
if lr
    tau_L = zeros(P*Nc*2,1); 
    for k = 1:P
        tau_k = [tau_stokes_x(Nc*(k-1)+1:Nc*k); tau_stokes_y(Nc*(k-1)+1:Nc*k)];
        tau_k = tau_k-L{1}*tau_k;
        tau_L(Nc*(k-1)+1:k*Nc)= tau_k(1:end/2);
        tau_L(Nc*(k-1)+P*Nc+1:k*Nc+P*Nc) = tau_k(end/2+1:end);
    end

    u_rhs = getVelocityField(rin,rcheck_b,tau_L(1:end/2),tau_L(end/2+1:end),[],[],[],[],[],[], []);
else
    u_rhs = matvec_mob_1B(mu_gmres,rin,rout,rcheck_b,rimage,nimage,q,UU,Y,L,pair_points,s,0,project_proxy);
end
S_0 = getRecompletionFlow(rin,rcheck_b,q,F,T); 
u_rhs = u_rhs-S_0; 
% We should get the same thing by just determining the TOTAL source
% strenghts (without the split) and evaluating directly.


disp('Surface residual')
diff_vec = u_rhs-u_lhs';
max_abs = max(abs(u_rhs(1:end/2)+1i*u_rhs(end/2+1:end)));
res = abs(diff_vec(1:end/2)+1i*diff_vec(end/2+1:end));
abs_res = max(res);
if max_abs > 0
    rel_vec = res/max_abs;
else
    rel_vec = res;
end
rel_res = max(rel_vec);
fprintf('Relative boundary error: %.3e\n', rel_res);
fprintf('Absolute boundary error: %.3e\n', abs_res);


if visualise
    if strcmp(surface_error_mode,'rel')
        boundary_err_plot = rel_vec;
        boundary_err_label = 'relative';
    else
        boundary_err_plot = res;
        boundary_err_label = 'absolute';
    end

    %% Visualise residual
    figure()
    subplot(1,2,1)
    plot(u_rhs);
    hold on
    plot(u_lhs);
    axis tight
    legend('RHS flow field','LHS flow field','interpreter','latex')
    subplot(1,2,2)
    semilogy(boundary_err_plot);
    legend(sprintf('%s boundary error', boundary_err_label),'interpreter','none')
    axis tight
    title(sprintf('%s: Boundary %s error', solver_name, boundary_err_label), 'interpreter','none')
    sgtitle([solver_name ': check solution flow field'],'interpreter','none')
    
    %% Visualise density
    lambda_x = lambdahat(1:length(rin));
    lambda_y = lambdahat(length(rin)+1:2*length(rin));

    lambda_image = lambdahat(2*length(rin)+1:end);

    figure()
    lambda_tot = vecnorm([lambda_y lambda_x],2,2);
    for k = 1:P
        scatter(real(rin((k-1)*Nc+1:k*Nc)),imag(rin((k-1)*Nc+1:k*Nc)),...
            20,log10(abs(lambda_tot((k-1)*Nc+1:k*Nc))),'filled');
        hold on
    end
    %Visualise extra singularities at image points
    if size(rimage,1)
        lambda_image_1 = lambda_image(1:end/2);
        lambda_image_2 = lambda_image(end/2+1:end);
        lambda_image_x1 = lambda_image_1(1:end/2);
        lambda_image_y1 = lambda_image_1(end/2+1:end);
        lambda_image_x2 = lambda_image_2(1:end/2);
        lambda_image_y2 = lambda_image_2(end/2+1:end);
        lambda_tot1 = vecnorm([lambda_image_x1 lambda_image_y1],2,2);
        lambda_tot2 = vecnorm([lambda_image_x2 lambda_image_y2],2,2);
       
        %scatter(real(rimage),imag(rimage),30,log10(abs(lambda_tot1)),'filled');
        %scatter(real(rimage-0.01i),imag(rimage-0.01i),30,log10(abs(lambda_tot2)),'filled');
%         else
%             scatter(real(0.9*(rimage-center)+center),imag(0.9*(rimage-center)+center),...
%             30,log10(abs(lambda_tot)),'filled');
%         end
    

    end
    c = colorbar;
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    hold on
    plot(real(rout),imag(rout),'k.')
    ylabel(c,'Magnitude of MFS coefficients. log10','interpreter','latex')
    set(c,'TickLabelInterpreter','latex')
    set(c,'TickLabelInterpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    axis off
    set(gcf,'color','w');
    title('Source solution 1-body','interpreter','latex')
    
end

end

function test_solve_mob
close all; 
images = 0; %images not needed for well separated particles
q = [0; 2.01]; %center coordinates
F = [1 0; 0 0]; %forces on the particles
T = [1; 1]; %torques on the particles
rads = [1; 1]; 
visualise = 1; 
lr = 0; %long range preconditioning
[UW,lambda_mob,it,gmres_tol,err] = solve_mob_1B(q,F,T,rads,images, lr, visualise);

%compare to a solution with image enhancement
images = 1; 
[UW_im,lambda_mob,it,gmres_tol,err_im] = solve_mob_1B(q,F,T,rads,images, lr, visualise);

str = sprintf('Relative residual with image enhancement: %1.2e vs without: %1.2e',err_im,err);
disp(str)

alignfigs(4);


end
