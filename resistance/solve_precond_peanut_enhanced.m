function [FT,lambda_proxy,it,gmres_tol,maxres] = solve_precond_peanut_enhanced(q,U,W,rads,delta_pair,N_peanut,visualise,gmres_tol,debug)
%SOLVE_PRECOND_PEANUT_ENHANCED Solves a 2D Stokes resistance problem with circular
%particles using a 2-body preconditioned MFS formulation. A
%fine grid enhanced with shielding Stokeslets near image points is used locally for every
%close pair of particles to resolve challenging interactions. This fine grid is
%compressed so that only coarse sources and collocation points are needed to 
%solve the problem iteratively, effectively preconditioning the system.
%
% Syntax:
%   [FT, lambda_proxy, it, gmres_tol, maxres] = solve_precond_peanut_enhanced(q, U, W, rads, delta_pair, N_peanut, visualise)
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
%   solve_2D_res              - 1-body preconditioned resistance solver
%   solve_res_precond_enhanced - 2-body preconditioner without peanut compression
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
solver_name = 'solve_precond_peanut_enhanced';

% Grid params
%Set coarse and fine grid. 

%Play with N_c, N_f, a (a_f). 
N_c = 60;  %100 better here? 
%N_c = 50;

N_f = 60;
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

opt.cmap = 0;
opt.P = P; 
opt.use_fmm_velocity = true; % set false to evaluate Stokeslet part with stokesletDirect

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
[UB_all,YB_all,UC_all,YC_all,Cmap,~] = getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt);

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
% hold on
% semilogy(resvec2); 

% Prepare for evaluating flow field in rcheck_b.
geom_eval = geom;
geom_eval.rcheck = rcheck_b;
[tau_stokes_x, tau_self_x, tau_beta_x,tau_cf_x,tau_stokes_y,tau_self_y,tau_beta_y,tau_cf_y,u_corr,rimage_k] = ...
    transform_peanut_stokes(tau,geom_eval,basis);


lambda_proxy = [tau_stokes_x; tau_stokes_y];

warning('Check what sources to report')

%Not yet adopted to random evaluation points. 

%% Do the evaluation of the flow in check points 
ftest_b = getVelocityField(rvec_in_c, rcheck_b, tau_stokes_x, tau_stokes_y);
rcheck_dom = [100+100i; -50+50i];
ftest = getVelocityField(rvec_in_c, rcheck_dom, tau_stokes_x, tau_stokes_y);


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

maxres = max(sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./max(sqrt(fb_x.^2+fb_y.^2))

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

% Then, due to all pair sources (fine-body + fine-image).
has_neigh = sort(unique(pairs(:)));
for i = 1:length(has_neigh)
    k = has_neigh(i);
    src_k = [rbase_in_f+q(k); rimage_k{k}];
    Kpair = getKmat2D(src_k,q(k));
    FT((k-1)*3+1:3*k) = FT((k-1)*3+1:3*k)+Kpair'*[tau_beta_x{k}; tau_beta_y{k}];
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
        Nother = singleLayer(rbase_in_c+q(i),rbase_out_f+q(p2),mu);
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
        Nother = singleLayer(rbase_in_c+q(p2),rbase_out_f+q(i),mu);
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

        Npair = singleLayer(rpair_fine,rcheck,mu);
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
        Npair = singleLayer(rpair_fine,[rbase_out_f+q(i); rbase_out_f+q(p2)],mu);
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
test = 1; 
delta_pair = 0.2; 
N_peanut = 400; 

if test == 1
    delta = 0.01;
    q = [0; 2+delta; (2+delta)*1i]; %center coordinates
    
    %or, instead, three circles in triangle - this DOES NOT converge to a good
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

    q = [0; 2+delta];
   % q = q+5;
    U = U(1:2,:);
    W = W(1:2); 
    rads = rads(1:2); 

    gmres_tol = 1e-10;
    debug = 0; 
    
    [FT1,lambda,it1,gmres_res, err1] = solve_precond_peanut_enhanced(q,U,W,rads,delta_pair,N_peanut,visualise,gmres_tol,debug);
    [FT2,lambda2,it2,gmres_res2, err2] = solve_res_precond_peanut(q,U,W,rads,delta_pair,N_peanut,visualise,0,gmres_tol,debug);
    [FT3,lambda,it3,gmres_res, err3] = solve_res_precond_enhanced(q,U,W,rads,delta_pair,0,visualise,gmres_tol,debug);
    
    F = [FT2(1:3:end) FT2(2:3:end)];
    T = FT2(3:3:end); 
    [UW,lambdahat,it1,gmres_mob, err_mob] = solve_mob_precond_peanut(q,F,T,rads,delta_pair,N_peanut,visualise);
    Ures = [U W]';
    
    
    % Compare to solution with pair corrections but without peanut compression
    visualise = 0; 
    [FTp,lambda,it2,gmres_tol, err2] = solve_res_precond_images(q,U,W,rads,delta_pair,0,visualise,gmres_tol,debug);
    
    str = sprintf('Two way error is %1.3e',norm(Ures(:)-UW));
    disp(str); 
    str = sprintf('Relative residual with 2-body precond: %1.2e vs with peanut compression: %1.2e\n Converging in %u resp % u iterations',err2,err1,it2,it1);
    disp(str)
    
    alignfigs(4);
else

    rng(9);
    P = 10;
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
    [FT1,lambda1,it1,gmres_tol,err1] = solve_precond_peanut_enhanced(q,U,W,rads,delta_pair,N_peanut,visualise);
    [FT2,lambda2,it2,gmres_tol,err2] = solve_res_precond_images(q,U,W,rads,delta_pair,visualise);
    
    str = sprintf('Relative residual with peanut compression: %1.2e vs pair preconditioner: %1.2e\\n Converging in %u resp %u iterations',err1,err2,it1,it2);
    disp(str)
    
    alignfigs;




end

end
