function [FT,lambda_proxy,it,gmres_tol,maxres] = solve_res_peanut_images(q,U,W,rads,delta_pair,N_peanut,visualise,lr,gmres_tol,debug,gmres_verbose)
%SOLVE_RES_PEANUT_IMAGES Solves a 2D Stokes resistance problem with circular
%particles using a 2-body preconditioned MFS formulation. A
%fine grid enhanced with approximate image points is used locally for every
%close pair of particles to resolve challenging interactions. This fine grid is
%compressed so that only coarse sources and collocation points are needed to 
%solve the problem iteratively, effectively preconditioning the system.
%
% Syntax:
%   [FT, lambda_proxy, it, gmres_tol, maxres] = solve_res_peanut_images(q, U, W, rads, delta_pair, N_peanut, visualise,lr)
%
% Inputs:
%   q          - Vector of length P, complex-valued center coordinates for the particles
%   U          - Px2 matrix of net force vectors (columns: x and y components)
%   W          - Px1 column vector of torques acting on the particles
%   rads       - Px1 vector of particle radii
%   delta_pair - Scalar threshold used to determine which particle pairs are considered close. For such pairs, a fine BVP is solved locally (a pair correction is built).
%   N_peanut   - Number of points on the peanut separation surface between
%                every close pair of particles. The peanut boundary is used
%                to map fine sources to effective coarse sources, giving the
%                same flow field exterior to the close pair of particles.
%   visualise  - Logical flag: plot the configuration and solution details
%   lr         - Flag for long-range preconditioning
%   gmres_tol  - Optional GMRES tolerance (default 1e-10)
%   debug      - Optional logical flag: build/draw dense matrix CC and its
%                eigenvalues for diagnostics (default false)
%
% Outputs:
%   FT         - 3P×1 vector of computed net forces and torques 
%   lambda_proxy     - Solution vector of coarse source strengths
%   it         - Number of GMRES iterations required
%   gmres_tol  - Set GMRES tolerance
%   maxres     - Maximum relative residual in a test (non-collocation) set of boundary nodes
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
%   solve_res_2B_images - 2-body preconditioner without peanut compression
%   solve_mob_peanut_images   - 2-body preconditioned mobility solver
%   with peanut compression
%
% To test: Call without arguments.
%
% Anna Broms, April 10, 2025

if nargin==0, test_solve_res; 
    return; end

if nargin < 9 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 10 || isempty(debug), debug = false; end
if nargin < 11 || isempty(gmres_verbose), gmres_verbose = 0; end

P = length(q); % number of particles

%% Checks 

assert(size(W,1)==P,'Wrong size of angular velocity vector')
assert(size(U,1)==P,'Wrong size of trans vel vector')
assert(size(U,2)==2,'Wrong size of trans vel vector, should contain x y coordinates')


%% SET PARAMS
%GMRES params
maxit = 800; 

if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

% Grid params
%Set coarse and fine grid. 

%Play with N_c, N_f, a (a_f). 
N_c = 80;  %100 better here? 
%N_c = 50;

N_f = 150;
%N_f = N_c; %debug

%N_c = 250; 
a_c = 1.2; %upsampling for coarse grid

%a_c = 5;
%a_c = 1.2; %results in horrible fourier decay if we 
           % evaluate on a particle close by in coarse colloc points 

a_f = 1.2; %upsampling factor for the fine grid

tol_c = 1e-10; %I think this works reasonably
tol_c = 1e-12;
%tol_c = 1e-10; %Curve moves closer to the surface -> smaller coeff 
%tol_c = 1e-16; %Curve moves further from surface -> larger coeff. 


s = [0 0 1 1 0 0 0]; %set type of singularities at image points [S R T D]
%s = [1 0 1 1]; %Other singularities? Currently not supported! But code can
%be changed!

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rp_c = max([1-sep_c,0.01]); %radius of proxy surface for coarse grid
Rp_f = max([1-sep_f,0.01]);  % and fine grid
%Rp_c = 0.8;
%Rp_f = Rp_c; %debug


%accumulation point, given Rp and delta. Closed formula from fixed point of reflection formula
accstop = (1-Rp_c)^2/Rp_c;  

if nargin < 6
    delta_pair = accstop; %We want to use the pair correction for all gaps smaller than delta_pair. (or accstop).
end

opt = get2Dparams(); 
opt.gmres_verbose = gmres_verbose;
opt.s = s;
opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c; 
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 
opt.rads = rads; 
opt.N_peanut = N_peanut;
opt.s = s; 
opt.mask = 0; %cutoff in long range preconditioning
opt.precomp = 1; %faster if evaluation of one body basis on fine grid is compted only once. 
% %Less storage required.
opt.pc = 1; %prepare grid to do pair corrections
opt.delta_pair = delta_pair; 

opt.cmap = 1;
opt.lr = lr; %long-range preconditioning 
opt.P = P; 

%% CREATE GRID
%Outer basic grid
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';
tout_f = linspace(0,2*pi,ceil(a_f*N_f)+1);
tout_f = tout_f(1:end-1)';
rbase_out_f = cos(tout_f)+1i*sin(tout_f);
rbase_out_c = cos(tout_c)+1i*sin(tout_c);
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rp_c*cos(tin)+Rp_c*1i*sin(tin); 
rvec_in_c = [];
for k = 1:P
    rvec_in_c = [rvec_in_c; q(k)+rbase_in_c];
end

%Construct image grid
%Will return only the basic outer grid, else refined outer grid 
%[rout, weights, rin, rimage, nimage, pair_points, pairs, rimage_pairs, refine, rin_base]
[rout,~,~,~,~,pair_points,pairs,rimage_vec,refine,rbase_in_f] = get2DImageGrid(q,rads,opt);


%Get pair basis
[UB_all,YB_all,UC_all,YC_all,Cmap,~,nimage] = getPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,[],[],[]);

%Get one-body pseduo inverse blocks -- enough to do this for single body.
[UU,YY] = getSelfPseudo(1,rbase_in_c,rbase_out_c);


%Visualise 1-body and pair-basis
%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,UU,YY,[],[],N_c, N_f,a_c,a_f,rads)

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


%% Construct check boundaries
%rcheck_dom = 100+100i; %some point far away

% Create new grid points, for which the accuracy of the solution is
% to be evaluated. 
rcheck_b = [];
n_bound = 803;
t = linspace(0,2*pi,n_bound)';
for k = 1:P
    rcheck_b = [rcheck_b; q(k)+rads(k)*(cos(t)+1i*sin(t))];
end

%% Experiment with left preconditioner based on deflation

if lr
    rin_c = []; 
    for k = 1:P
        rin_c = [rin_c; rbase_in_c+q(k)];
    end
    [Sinv,Zi,Yi,db] = get_long_range_precond_peanut(q,rin_c,rout,rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,opt);

    %[Sinv,Zi,Yi,db] = get_long_range_precond_mu(q,rin_c,rout,opt);
    opt.db = db;
    mu_coarse = getCoarseMu(fout,Sinv,Zi,Yi,db,P,opt.N_c,opt.a_c);

    %Try to evaluate: 
    res1 = matvec_res_peanut(mu_coarse,rbase_in_c,rbase_in_f,rvec_in_c,rbase_out_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,0);
    
    %pair_points = [zeros(P,1) opt.a_c*opt.N_c*ones(P,1)];
    %res2 = matvec_res_Stokes(mu_coarse,rvec_in_c,rout,[],[],q,UU,YY,pair_points,s);
end

%% Solve system

% Build the matrix to inspect conditioning/eigenvalues if requested.
if debug
    matvec_debug = 0;
    x = zeros(2*length(rout),1);
    tic
    ncols = 2*length(rout);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        uu = matvec_res_peanut(x,rbase_in_c,rbase_in_f,rvec_in_c,rbase_out_f,...
            refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,matvec_debug);
        CC(:,k) = uu;
    end
    toc
    figure(14);
    clf; 
    imagesc(log10(abs(CC)))
    colorbar
    title([solver_name ': log_{10} |CC|'],'interpreter','none')
    cc = skeel(CC);
    fprintf('Estimated condition number of system matrix: %1.3e \n',cc);
    figure();
    [~,D] = eig(CC);
    D = diag(D); 
    plot(real(D),imag(D),'b+')
    title([solver_name ': eigenvalues of CC'],'interpreter','none')
end

% To test with Krylov precond, do something like
%[tau1, e1, precond] = precond_gmres(@(x) matvec_res_peanut(x,rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,opt,rvec_out,q,UU,YY,pairs,UB_all,YB_all,DC_all, YC_all,debug), fout, zeros(2*size(rvec_out,1),1), 2*size(rvec_out,1), gmres_tol, precond);
%fprintf("iterations = %d\n", length(e1))

disp(' == Solving... == ');
if lr
   Pf = applyPmat_peanut(fout,rin_c,rout,Sinv,q,Zi,Yi,rbase_in_c,...
        rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,opt,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap);
   %Pf = applyPmat_peanut(fout,rin_c,rout,Sinv,Nx,Ny,Mx,Zi,Yi,opt);    
   [tau,it,resvec,real_res] = helsing_gmres(@(x) lr_matvec_2D_peanut(x,rin_c,rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Sinv,Zi,Yi),Pf,2*size(rout,1),maxit,gmres_tol,opt,rout);   
else                                                                                 
   matvec_debug = 0;
   [tau,it,resvec,real_res] = helsing_gmres(@(x) matvec_res_peanut(x,rbase_in_c,rbase_in_f,rvec_in_c,rbase_out_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,matvec_debug),fout,2*size(rout,1),maxit,gmres_tol,opt,rout);
end

%[tau,flag,relres,iter,resvec2] = gmres(@(x) matvec_2D_pairprecond3(x,rbase_in_c,rbase_in_f,rbase_out_f,rvec_out,q,UU,YY,B,pairs,A,Uf,Yf,Ncf,Upf,Ypf),fout,[],gmres_tol,maxit);
figure()
semilogy(resvec); 
title('GMRES convergence with peanut compression, resistance', 'interpreter','latex')

if visualise
    %check residual
    matvec_debug = 0;
    restot = (matvec_res_peanut(tau,rbase_in_c,rbase_in_f,rvec_in_c,rbase_out_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,matvec_debug)-fout);
    figure()
    semilogy(abs(restot))
    title('Res at colloc points, peanut resistance')
end
% hold on
% semilogy(resvec2); 



disp(' == Postprocessing == ');
if lr

    mu_mapped = applyQmat_peanut(tau,rvec_in_c,rout,Sinv,Zi,Yi,opt,...
        rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap);

    mu_tot = mu_mapped+mu_coarse; 

    [tau_stokes_x, tau_self_x, tau_beta_x,tau_stokes_y,tau_self_y,tau_beta_y,u_corr] = transform_peanut(mu_tot,...
        rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,opt,rout,rcheck_b,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,0,0);


    N_large = N_c*a_c;
    PM = length(rout);
    tau_coarse_x = zeros(N_c*P,1);
    tau_coarse_y = zeros(N_c*P,1);
    for i = 1:P
        coarse_ind = (i-1)*N_c+1:i*N_c; 

        %Retrieve self evaluation blocks    
        tau_particle_x = mu_coarse((i-1)*N_large+1:N_large*i);
        tau_particle_y = mu_coarse(PM+(i-1)*N_large+1:PM+N_large*i);

        %U{i}, Y{i}
        step1 = UU{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other?
        tau_mapped = YY{1}*step1; %this is the mapped density for this particle to throw in to the kernel
        tau_coarse_x(coarse_ind) = tau_mapped(1:N_c);
        tau_coarse_y(coarse_ind) = tau_mapped(N_c+1:end);

    end
    % 
    % tau_stokes_x2 = tau_stokes_x; 
    % tau_stokes_y2 = tau_stokes_y; 
    % 
    % %think about the modification to make here
    % %tau_stokes = applyQmat([tau_stokes_x; tau_stokes_y],rvec_in_c,rout,Sinv,[],[],[],Zi,Yi,opt);
    % tau_stokes = applyQmat_peanut(tau,[tau_stokes_x; tau_stokes_y],rvec_in_c,rout,Sinv,Zi,Yi,opt,...
    %     rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap);
    % tau_stokes_x = tau_stokes(1:end/2);%tau_coarse(1:end/2);
    % tau_stokes_y = tau_stokes(end/2+1:end);%tau_coarse(end/2+1:end);
    % 
    % %Now, what if I map this back to mu?
    % % need single body forward SLP
    % S = stokSLPmat(rbase_in_c,rout(1:N_c*a_c),1);
    % tau_new = zeros(size(tau));
    % for i = 1:P
    % 
    %     mu_i = S*[tau_stokes_x((i-1)*N_c+1:i*N_c); tau_stokes_y((i-1)*N_c+1:i*N_c)];
    %     tau_new((i-1)*N_c*a_c+1:N_c*a_c*i) = mu_i(1:end/2); 
    %     tau_new((i-1)*N_c*a_c+N_c*a_c*P+1:N_c*a_c*i+N_c*a_c*P) = mu_i(end/2+1:end); 
    % 
    %     %Only modify self-blocks:
    % 
    % 
    % end
    % %Prepare for evaluating flow field and evaluate in new points rcheck_dom and rcheck_b
    % [tau_stokes_x, tau_self_x, tau_beta_x,tau_stokes_y,tau_self_y,tau_beta_y,u_corr] = transform_peanut(tau_new,...
    %     rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,opt,rout,rcheck_b,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,0,0);
       
else
    %Prepare for evaluating flow field and evaluate in new points rcheck_dom and rcheck_b
    [tau_stokes_x, tau_self_x, tau_beta_x,tau_stokes_y,tau_self_y,tau_beta_y,u_corr] = transform_peanut(tau,...
        rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,opt,rout,rcheck_b,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,0,0);
end

% if lr 
%     tau_stokes_x = tau_stokes_x+tau_coarse_x;
%     tau_stokes_y = tau_stokes_y+tau_coarse_y;
% 
%     tau_self_x = tau_self_x+tau_coarse_x;
%     tau_self_y = tau_self_y+tau_coarse_y;
% end


lambda_proxy = [tau_stokes_x; tau_stokes_y];

warning('Check what sources to report')

%Not yet adopted to random evaluation points. 

%% Do the evaluation of the flow in check points 
ftest_b = getVelocityField(rvec_in_c, rcheck_b, tau_stokes_x, tau_stokes_y, [], [], ...
                           [],[],[],[],[]);
rcheck_dom = [100+100i; -50+50i];
ftest = getVelocityField(rvec_in_c, rcheck_dom, tau_stokes_x, tau_stokes_y,[], [], ...
                           [], [], [],[], [])


ftest_b = ftest_b+u_corr; 

%Compute error in the value at the boundary
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

maxres = max(sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./max(sqrt(fb_x.^2+fb_y.^2));
fprintf('Max surf rel res at new nodes %.3e\n', maxres);

% % For debuggin: 
% fbound_x2 = ftest_b2(1:length(rcheck_b));
% fbound_y2 = ftest_b2(length(rcheck_b)+1:end);
% err2 = max(sqrt((fb_x+fbound_x2).^2+(fb_y+fbound_y2).^2))./max(sqrt(fb_x.^2+fb_y.^2))

%% Determine force and torque
%from coarse sources
K = getKmat2D(rbase_in_c,0);
FT = zeros(3*P,1); 
for k= 1:P
    FT((k-1)*3+1:3*k) = K'*[tau_self_x((k-1)*N_c+1:k*N_c); tau_self_y((k-1)*N_c+1:k*N_c)];
end

% %Then, due to all fine sources.
Kf = getKmat2D(rbase_in_f,0);
has_neigh = sort(unique(pairs(:)));
for i = 1:length(has_neigh)
    k = has_neigh(i); 
    FT((k-1)*3+1:3*k) = FT((k-1)*3+1:3*k)+Kf'*[tau_beta_x((k-1)*N_f+1:k*N_f); 
        tau_beta_y((k-1)*N_f+1:k*N_f)];
end

% % Alternative way to account for all the fine sources
% for i = 1:length(has_neigh)
%     k = has_neigh(i); 
%     FT2((k-1)*3+1:3*k) = K'*[tau_stokes_x((k-1)*N_c+1:k*N_c); 
%         tau_stokes_y((k-1)*N_c+1:k*N_c)];
% end
% 
% %extract force and torque separately
% for k = 1:P
%     F(k) = FT((k-1)*3+1) + 1i*FT((k-1)*3+2);
%     T(k) = FT(3*k); 
%     F2(k) = FT2((k-1)*3+1) + 1i*FT2((k-1)*3+2);
%     T2(k) = FT2(3*k); 
% end

if visualise

    figure()
    plot(ftest_b)
    hold on
    plot([fb_x; fb_y])
    title('Velocity rhs and lhs, peanut resistance')

    %Check each component
%     figure(9)
%     subplot(2,2,1)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fb_x-fbound_x)),30,log10(abs(fb_x-fbound_x)),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('error in x velocity')
% 
% 
%     subplot(2,2,2)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fb_y-fbound_y)),30,log10(abs(fb_y-fbound_y)),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('error in y velocity')
% 
%     % Visualise the actual velocity
%     subplot(2,2,3)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fb_x)),30,log10(abs(fb_x)),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('x velocity MFS')
% 
%     subplot(2,2,4)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fbound_x)),30,log10(abs(fbound_x)),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('x velocity reference')
%
%     sgtitle('Error on boundary','interpreter','latex')

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

    %% Visualise density
    figure()
    subplot(1,2,1)
    semilogy(abs(lambda_proxy))
    hold on
    subplot(1,2,2)
    plot(lambda_proxy)
    sgtitle('Source strengths $\lambda$ with peanut compression, resistance', 'interpreter','latex')

    

    % lambda_x = lambda(1:length(rin_f));
    % lambda_y = lambda(length(rin_f)+1:2*length(rin_f));
    % 
    % figure(12)
    % p = length(rin_f)/length(q);
    % 
    % clf; 
    % lambda_tot = vecnorm([lambda_y lambda_x],2,2);
    % for k = 1:length(q)
    %     scatter(real(rin_f((k-1)*p+1:k*p)),imag(rin_f((k-1)*p+1:k*p)),...
    %         30,log10(abs(lambda_tot((k-1)*p+1:k*p))),'filled');
    %     hold on
    % end
    % %Visualise extra singularities at image points
    % 
    % c = colorbar;
    % axis equal
    % view(0,90)
    % grid off
    % set(gca,'xtick',[])
    % set(gca,'ytick',[])
    % hold on
    % plot(real(rout),imag(rout),'k.')
    % ylabel(c,'Magnitude of MFS coefficients. log10','interpreter','latex')
    % set(c,'TickLabelInterpreter','latex')
    % set(gca,'TickLabelInterpreter','latex')
    % axis off
    % set(gcf,'color','w');
end

end


%OLD FUNCTIONS... 

function doPairBasisTest(Upf,Ypf,i,p2,q,U,Y,rbase_in_c,rbase_out_f,rpair_fine)
%Test pair basis with a smooth coarse density
        mu = 1; 

        %generate density
        t = linspace(0,2*pi,ceil(1.2*size(rbase_in_c,1))+1);
        t = t(1:end-1)';
        tau_1 = sum([exp(1i*t),exp(2*1i*t),exp(3*1i*t),exp(4*1i*t),exp(4*1i*t)]*rand(5),2);
        tau_2 = sum([exp(1i*t),exp(2*1i*t),exp(3*1i*t),exp(4*1i*t),exp(4*1i*t)]*rand(5),2);

        tau_p1_x = real(tau_1);
        tau_p1_y = imag(tau_1); 
        tau_p2_x = real(tau_2);
        tau_p2_y = imag(tau_2);

        figure()
        plot(tau_p1_x);
        hold on
        plot(tau_p1_y);
        plot(tau_p2_x);
        plot(tau_p2_y);

        %Get A1:
        step1 = U'*[tau_p1_x;tau_p1_y]; %here I assume x and y follow each other?
        tau_mapped = Y*step1;
        Nother = stokSLPmat(rbase_in_c+q(i),rbase_out_f+q(p2),mu);
        R1 = -Nother*tau_mapped; %not correct. Should read off on particle 2
        block = R1(1:end/2);
    
        A2 = [zeros(size(block)); block; zeros(size(block)); R1(end/2+1:end)]; 
        %A1 = [tau_particle_x; zeros(size(tau_particle_x)); tau_particle_y; zeros(size(tau_particle_x))];
        pair_mapped = Upf{i,p2}'*A2;

        %pair_mapped2 = Upf{i,p2}'*A2;
        tau_mapped = Ypf{i,p2}*pair_mapped;

        %Get A2
        step1 = U{p2}'*[tau_p2_x;tau_p2_y]; %here I assume x and y follow each other?
        mapped = Y{p2}*step1;
        Nother = stokSLPmat(rbase_in_c+q(p2),rbase_out_f+q(i),mu);
        R1 = -Nother*mapped;
        block = R1(1:end/2);
        A1 = [R1(1:end/2); zeros(size(block)); R1(end/2+1:end); zeros(size(block))];
        %A2 = [zeros(size(block)); block; zeros(size(block)); R1(end/2+1:end)];            
        %A2 = [R1(1:end/2); zeros(size(block)); R1(end/2+1:end); zeros(size(block))];
        pair_mapped = Upf{i,p2}'*A1;
        tau_mapped2 = Ypf{i,p2}*pair_mapped;

        tau_tot = tau_mapped+tau_mapped2; 

        %Now apply. Want to compute what this is exterior to the two
        %particles in the pair
        x = linspace(-2,5);
        y = linspace(-2,5);
        [X,Y] = meshgrid(x,y); 
        rcheck = X(:)+1i*Y(:);
        ind1 = find(abs(rcheck-q(i))<1);
        ind2 = find(abs(rcheck-q(p2))<1);
        ind_keep_1 = setdiff(1:size(rcheck,1),ind1);
        ind_keep_2 = setdiff(ind_keep_1,ind2);
        %rcheck = rcheck(ind_keep_2);

        Npair = stokSLPmat(rpair_fine,rcheck,mu);
        ucheck = Npair*tau_tot;
      %  ucheck(ind1) =  nan;
      %  ucheck(ind2) = nan; 
        ucheck = ucheck(1:end/2)+1i*ucheck(end/2+1:end);
        ucheck(ind1) = nan;
        ucheck(ind2) = nan;


        figure()
        surfir(real(rcheck),imag(rcheck),abs(ucheck));
        colorbar
        axis equal
        hold on
        plot(real(q),imag(q),'k+')
        view(0,90);

        figure()
        surfir(real(rcheck),imag(rcheck),log10(abs(ucheck)));
        colorbar
        axis equal
        hold on
        plot(real(q),imag(q),'k+')
        view(0,90);

        %Also, compute residuals from the solve steps above.
        %ucheck is the velocity computed with the pair basis?
        Npair = stokSLPmat(rpair_fine,[rbase_out_f+q(i); rbase_out_f+q(p2)],mu);
        ulhs = Npair*tau_mapped;
        urhs = A2;

        ulhs2 = Npair*tau_mapped2;
        urhs2 = A1;

        norm(ulhs-A2,inf) %seems pretty accurate... 
        norm(ulhs2-A1,inf)

%         figure(10)
%         clf;
%         plot((A1(1:end/4)-ulhs2(1:end/4))./A1(1:end/4))
%         hold on
%         plot((A1(end/2+1:3*end/4)-ulhs2(end/2+1:3*end/4))./A1(end/2+1:3*end/4))
%         plot((A2(end/4+1:end/2)-ulhs(end/4+1:end/2))./A2(end/4+1:end/2))
%         plot((A2(3*end/4+1:end)-ulhs(3*end/4+1:end))./A2(3*end/4+1:end))

        figure(10)
        clf;
        plot((A1(1:end/4)-ulhs2(1:end/4)))
        hold on
        plot((A1(end/2+1:3*end/4)-ulhs2(end/2+1:3*end/4)))
        plot((A2(end/4+1:end/2)-ulhs(end/4+1:end/2)))
        plot((A2(3*end/4+1:end)-ulhs(3*end/4+1:end)))


end

function test_solve_res

close all; 
test = 2; 
delta_pair = 0.2; 
N_peanut = 400; 

if test == 1
    delta = 0.001;
    q = [0; 2+delta; (2+delta)*1i]; %center coordinates
    
    %or, instead, three cireles in triangle - this DOES NOT converge to a good
    %sol! 
    % x = 1+delta/2;
    % y = sqrt((2+delta)^2-(1+delta/2)^2);
    % q = [0; 2+delta; x+1i*y];
    
    
    U = [1 0; 0 0; 0 1]; %translational velocities for particles
    W = [1; 1; 1]; %angular velocities
    % U = U*1e-5;
    % W = W*1e-5;
    rads = [1; 1; 1]; 
    visualise = 1; 
    
    [FT,lambda,it1,gmres_res, err1] = solve_res_peanut_images(q,U,W,rads,delta_pair,N_peanut,visualise);
    F = [FT(1:3:end) FT(2:3:end)];
    T = FT(3:3:end); 
    [UW,lambdahat,it1,gmres_mob, err_mob] = solve_mob_peanut_images(q,F,T,rads,delta_pair,N_peanut,visualise);
    Ures = [U W]';
    
    
    % Compare to solution with pair corrections but without peanut compression
    visualise = 0; 
    [FTp,lambda,it2,gmres_tol, err2] = solve_res_2B_images(q,U,W,rads,delta_pair,visualise);
    
    str = sprintf('Two way error is %1.3e',norm(Ures(:)-UW));
    disp(str); 
    str = sprintf('Relative residual with 2-body precond: %1.2e vs with peanut compression: %1.2e\n Converging in %u resp % u iterations',err2,err1,it2,it1);
    disp(str)
    
    alignfigs(4);
else

    rng(9);
    P = 2;
    delta = 0.001; %P = 5
    x = 1+delta/2;
    y = sqrt((2+delta)^2-(1+delta/2)^2);
    q = [0; 2+delta; x+1i*y];
   % P = length(q); 
    
    visualise = 1; 
    %delta = 1;
    q = grow_cluster(P,delta,2);
   % q = [q; -6+1.5i; -2-4i]; P = P+2;
    %q = q([1,2,4],:); P = 3; 
    U = rand(P,2); W = rand(P,1); rads = ones(P,1);
    %W = zeros(P,1); 
    lr = 20; 
   % lr = 0; 
    images = 0; 
   % [FT,lambda,it,gmres_tol,err] = solve_res_1B(q,U,W,rads,images, lr,visualise);
    [FT1,lambda1,it1,gmres_tol,err1] = solve_res_peanut_images(q,U,W,rads,delta_pair,N_peanut,visualise,lr);
    [FT2,lambda2,it2,gmres_tol,err2] = solve_res_peanut_images(q,U,W,rads,delta_pair,N_peanut,visualise,0);
    
    str = sprintf('Relative residual with 1-body precond: %1.2e vs 2-body: %1.2e\n Converging in %u resp % u iterations',err1,err2,it1,it2);
    disp(str)
    
    alignfigs;




end

end
