function [FT,lambda_proxy,it,gmres_tol,maxres] = solve_res_peanut_enhanced(q,U,W,rads,delta_pair,N_peanut,visualise,gmres_tol,debug)
%SOLVE_PRECOND_PEANUT_ENHANCED Solves a 2D Stokes resistance problem with circular
%particles using a 2-body preconditioned MFS formulation. A
%fine grid enhanced with shielding Stokeslets near image points is used locally for every
%close pair of particles to resolve challenging interactions. This fine grid is
%compressed so that only coarse sources and collocation points are needed to 
%solve the problem iteratively, effectively preconditioning the system.
%
% Syntax:
%   [FT, lambda_proxy, it, gmres_tol, maxres] = solve_res_peanut_enhanced(q, U, W, rads, delta_pair, N_peanut, visualise)
%
% Inputs:
%   q          - Vector of length P, complex-valued center coordinates for the particles
%   U          - Px2 matrix of translational velocities (columns: x and y components)
%   W          - Px1 column vector of angular velocities
%   rads       - Px1 vector of particle radii
%   delta_pair - Scalar threshold used to determine which particle pairs are considered close. For such pairs, a fine BVP is solved locally (a pair correction is built).
%   N_peanut   - Number of points on the peanut separation surface between
%                every close pair of particles. The peanut boundary is used
%                to map fine sources to effective coarse sources, giving the
%                same flow field exterior to the close pair of particles.
%   visualise  - Logical flag: plot the configuration and solution details
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
%   solve_res_2B_enhanced     - 2-body preconditioner without peanut compression
%   solve_mob_peanut_enhanced   - 2-body preconditioned mobility solver
%   with peanut compression
%
% To test: Call without arguments.
%
% Anna Broms, Feb 23, 2026

if nargin==0, test_solve_res; 
    return; end

if nargin < 8 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 9 || isempty(debug), debug = false; end

P = length(q); % number of particles

%% Checks 

assert(size(W,1)==P,'Wrong size of angular velocity vector')
assert(size(U,1)==P,'Wrong size of trans vel vector')
assert(size(U,2)==2,'Wrong size of trans vel vector, should contain x y coordinates')


%% SET PARAMS
%GMRES params
maxit = 800; 
solver_name = 'solve_res_peanut_enhanced';

% Grid params
%Set coarse and fine grid. 

%Play with N_c, N_f, a (a_f). 
N_c = 150;  %100 better here? 
%N_c = 50;

N_f = 150;
%N_f = N_c; %debug

%N_c = 250; 
a_c = 1.2; %upsampling for coarse grid

%a_c = 5;
%a_c = 1.2; %results in horrible fourier decay if we 
           % evaluate on a particle close by in coarse colloc points 

a_f = 1.2; %upsampling factor for the fine grid

tol_c = 1e-9; %I think this works reasonably
%tol_c = 1e-12;
%tol_c = 1e-10; %Curve moves closer to the surface -> smaller coeff 
%tol_c = 1e-16; %Curve moves further from surface -> larger coeff. 

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rp_c = max([1-sep_c,0.01]); %radius of proxy surface for coarse grid
Rp_f = max([1-sep_f,0.01]);  % and fine grid
%Rp_c = 0.8;
%Rp_f = Rp_c; %debug


%accumulation point, given Rp and delta. Closed formula from fixed point of reflection formula
accstop = (1-Rp_c)^2/Rp_c;  

if nargin < 5
    delta_pair = accstop; %We want to use the pair correction for all gaps smaller than delta_pair. (or accstop).
end

opt = get2Dparams(); 
opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c; 
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 
opt.Nclust = 200;
opt.N_peanut = N_peanut; 
opt.precomp = 1; %faster if evaluation of one body basis on fine grid is compted only once. 
% %Less storage required.
opt.pc = 1; %prepare grid to do pair corrections
opt.delta_pair = delta_pair; 

opt.cmap = 1;
opt.P = P; 
opt.use_fmm_velocity = true; % set false to evaluate Stokeslet part with stokesletDirect

opt.bndry_vel = 0; 
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
 
[~, ~, ~, rimage_vec, refine,pairs] = getEnhancedGrid(q, opt);


%Get pair basis                                               
opt.show_counter = true;
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU] = getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,[]);

%Get one-body pseduo inverse blocks -- enough to do this for single body.
[UU,YY] = getSelfPseudo(1,rbase_in_c,rbase_out_c);

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

%Visualise 1-body and pair-basis
%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,[],refine,UB_all,YB_all,UU,YY,[],[],N_c,N_f,a_c,a_f,rads)

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

%% Solve system

% Build the matrix to inspect conditioning/eigenvalues if requested.
if debug
    x = zeros(2*length(rout),1);
    tic
    for k = 1:2*length(rout)
        k
        x(:) = 0; 
        x(k) = 1; 
        uu = matvec_2D_peanut_enhanced(x,geom,basis);
        CC(:,k) = uu;
    end
    toc
    figure(14);
    clf; 
    imagesc(log10(abs(CC)))
    colorbar
    title([solver_name ': log_{10} |CC|'],'interpreter','none')
    skeel(CC)

    figure();
    [~,D] = eig(CC);
    D = diag(D); 
    plot(real(D),imag(D),'b+')
    title([solver_name ': eigenvalues of CC'],'interpreter','none')
end

[tau,it,resvec,real_res] = helsing_gmres( ...
    @(x) matvec_2D_peanut_enhanced(x,geom,basis), ...
    fout,2*size(rout,1),maxit,gmres_tol,1,rout);

figure()
semilogy(resvec); 
title('GMRES convergence with peanut compression, resistance', 'interpreter','latex')

if visualise
    %check residual
    restot = matvec_2D_peanut_enhanced(tau,geom,basis)-fout;
    figure()
    semilogy(abs(restot))
    title('Res at colloc points, peanut resistance')
end

%% Reconstruct sources. If the boundary velocities on each disc is sought, the fine sources are needed

if opt.bndry_vel
    % Create new grid points, for which the accuracy of the solution is
    % to be evaluated. 
    rcheck_b = [];
    n_bound = 803;
    t = linspace(0,2*pi,n_bound)';
    for k = 1:P
        rcheck_b = [rcheck_b; q(k)+rads(k)*(cos(t)+1i*sin(t))];
    end

    % Prepare for evaluating flow field in rcheck_b.
    geom_eval = geom;
    geom_eval.rcheck = rcheck_b;
    [lam_c_x, lam_self_x, lam_f_x,lam_c_y,lam_self_y,lam_f_y,u_corr,rimage_k] = ...
        transform_peanut_stokes(tau,geom_eval,basis);

    %% Do the evaluation of the flow in check points 
    ftest_b = getVelocityField(rvec_in_c, rcheck_b, lam_c_x, lam_c_y);
    ftest_b = ftest_b+u_corr; % Apply on pair corrections
    
    %% Compute error at these boundary nodes 
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
else
    geom_eval = geom;
    [lam_c_x, lam_self_x, ~,lam_c_y,lam_self_y,~,~,rimage_k] = ...
        transform_peanut_stokes(tau,geom_eval,basis);
    maxres = nan; 
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
if opt.cmap
    % Use precomputed map 
    for pair_it = 1:size(pairs,1)
        i = pairs(pair_it,1);
        p2 = pairs(pair_it,2);

        coarse_i = (i-1)*N_c+1:i*N_c;
        coarse_p2 = (p2-1)*N_c+1:p2*N_c;
        rhs_pair = [lam_self_x(coarse_i); lam_self_x(coarse_p2); ...
                    lam_self_y(coarse_i); lam_self_y(coarse_p2)];

        % Sign matches the fine-source postprocessing path via lam_f_{x,y}.
        FT_pair = Cmap_FU{i,p2}*rhs_pair; % [Fx_i; Fy_i; T_i; Fx_p2; Fy_p2; T_p2]
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


if visualise

    if opt.bndry_vel
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
    hold on
    subplot(1,2,2)
    plot(lambda_proxy)
    sgtitle('Source strengths $\lambda$ with peanut compression, resistance', 'interpreter','latex')

   
end

end


function test_solve_res

close all; 
test = 2; 
delta_pair = 0.2; 
N_peanut = 400; 

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
    rads = [1; 1; 1]; 
    visualise = 0; 

    q = [0; 2+delta];
    %q = q+5;
    U = U(1:2,:);
    W = W(1:2); 
    rads = rads(1:2); 

    gmres_tol = 1e-10;
    debug = 1; 
    
    [FT1,lambda,it1,gmres_res, err1] = solve_res_peanut_enhanced(q,U,W,rads,delta_pair,N_peanut,visualise,gmres_tol,debug);
    debug = 0; 
    [FT2,lambda2,it2,gmres_res2, err2] = solve_res_peanut_images(q,U,W,rads,delta_pair,N_peanut,visualise,0,gmres_tol,debug);
    [FT3,lambda,it3,gmres_res, err3] = solve_res_2B_enhanced(q,U,W,rads,delta_pair,0,visualise,gmres_tol,debug);
    
    F = [FT2(1:3:end) FT2(2:3:end)];
    T = FT2(3:3:end); 
    [UW,lambdahat,it1,gmres_mob, err_mob] = solve_mob_peanut_images(q,F,T,rads,delta_pair,N_peanut,visualise);
    Ures = [U W]';
    
    
    % Compare to solution with pair corrections but without peanut compression
    visualise = 0; 
    [FTp,lambda,it2,gmres_tol, err2] = solve_res_2B_images(q,U,W,rads,delta_pair,0,visualise,gmres_tol,debug);
    
    str = sprintf('Two way error is %1.3e',norm(Ures(:)-UW));
    disp(str); 
    str = sprintf('Relative residual with 2-body precond: %1.2e vs with peanut compression: %1.2e\n Converging in %u resp % u iterations',err2,err1,it2,it1);
    disp(str)
    
    alignfigs(4);
else

    rng(9);
    P = 4;
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
    gmres_tol = 1e-7; 
    debug = 0; 
    [FT1,lambda1,it1,gmres_tol,err1] = solve_res_peanut_enhanced(q,U,W,rads,delta_pair,N_peanut,visualise,gmres_tol,debug);
    [FT2,lambda2,it2,gmres_tol,err2] = solve_res_2B_images(q,U,W,rads,delta_pair,0,visualise);
    
    str = sprintf('Relative residual with peanut compression: %1.2e vs pair preconditioner: %1.2e\\n Converging in %u resp %u iterations',err1,err2,it1,it2);
    disp(str)
    
    alignfigs;




end

end
