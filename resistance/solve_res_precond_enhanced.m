function [FT,lambda,it,gmres_tol,maxres] = solve_res_precond_enhanced(q,U,W,rads,delta_pair,lr,visualise,gmres_tol,debug)
%SOLVE_RES_PRECOND_ENHANCED Solves a 2D Stokes resistance problem with circular
%particles using MFS with 2-body preconditioning. To resolve
% challenging close interactions, a fine 2-body BVP is solved for fine
% source strengths on each close pair of particles. The fine sources,
% including shielding of image points, correct the representation
% obtained from a coarse grid, effectively preconditioning the system.
%
% Syntax:
%   [FT, lambda, it, gmres_tol, maxres] = solve_res_precond_enhanced(q, U, W, rads, delta_pair, lr, visualise)
%
% Inputs:
%   q          - Vector of length P, complex-valued center coordinates for the particles
%   U         - Px2 matrix of translational velocities (columns: x and y components)
%   W         - Px1 column vector of angular velocities
%   rads       - Px1 vector of particle radii
%   delta_pair - Scalar threshold used to determine which particle pairs are considered close. For such pairs, a fine BVP is solved locally (a pair correction is built).
%   lr         - long range preconditioning 
%   visualise  - Logical flag: plot the configuration and solution details
%   gmres_tol  - Optional GMRES tolerance (default 1e-10)
%   debug      - Optional logical flag: build/draw dense matrix CC and its
%                eigenvalues for diagnostics (default false)
%
% Outputs:
%   FT         - 3P×1 vector of computed net forces and torques 
%   lambda     - Solution vector of source strengths
%   it         - Number of GMRES iterations required
%   gmres_tol  - Set GMRES tolerance
%   maxres     - Maximum relative residual in a test (non-collocation) set of boundary nodes
%
% Description:
%   This function applies a 2-body preconditioner (using pair corrections via local fine BVPs)
%   to improve solution accuracy for near-contact particle configurations.
%   The FMM is used for Stokeslet evaluation.
%
% Notes:
%   - This is an MFS generalisation of the Cheng-Greengard idea (1998/2000)
%   - Intended to match solve_mob_precond_images, but applied to a resistance problem rather than mobility.
%
% See also:
%   solve_2D_res              - 1-body preconditioned resistance solver
%   solve_res_precond_peanut - 2-body preconditioner with peanut
%                              compression (in this version of the algorithm, only coarse sources
%                              are needed at the solve stage).
%   solve_mob_precond_images   - 2-body preconditioned mobility solver
%
% To test: Call without arguments.
%
% Anna Broms, Feb 13, 2025

if nargin==0, test_solve_res; 
    return; end

if nargin < 8 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 9 || isempty(debug), debug = false; end

P = length(q);

%% Checks

assert(size(W,1)==P,'Wrong size of angular velocity vector')
assert(size(U,1)==P,'Wrong size of trans vel vector')
assert(size(U,2)==2,'Wrong size of trans vel vector, should contain x y coordinates')


%% SET PARAMS
%GMRES params
maxit = 800; 
solver_name = 'solve_res_precond_enhanced';

% Grid params
P = length(q); 

opt = get2Dparams();


%Play with N_c, N_f, a (a_f). 
N_c = 60;  %100 better here? 
%N_c = 80; 
%N_c = 150; 
%N_c = 100; 
N_f = 150; 
N_f = 60;
%N_f = N_c; %debug

%N_c = 250; 
a_c = 1.2;
 
%a = 2; 
a_f = 1.2; %upsampling factor for the fine grid 

% Determine separation to proxy surface.
tol_c = 1e-12; %I think this works reasonably

%tol_c = 1e-10; %Curve moves closer to the surface -> smaller coeff 
%tol_c = 1e-16; %Curve moves further from surface -> larger coeff. 

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rp_c = max([1-sep_c,0.01]); %radius of proxy surface for coarse grid
Rp_f = max([1-sep_f,0.01]);  % and fine grid
%Rp_f = Rp_c; %debugs

 
%accumulation point, given Rp and delta. Closed formula from fixed point of reflection formula
accstop = (1-Rp_c)^2/Rp_c;  

if nargin < 5
    delta_pair = accstop; %We want to use the pair correction for all gaps smaller than delta_pair. (or accstop).
end
 
opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c;
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 

opt.P = P; 
opt.N_peanut = 0; 
opt.precomp = 1; %faster if evaluation of one body basis on fine grid is compted only once. 
% %Less storage required.
opt.pc = 1; %prepare grid to do pair corrections
opt.delta_pair = delta_pair; 
opt.lr = lr; %using long range preconditioning? 
opt.Nclust = 200; %trial nodes at ellipse segments. 


%% CREATE GRID
%Outer basic grid
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';
rbase_out_c = rads(1)*cos(tout_c)+1i*rads(1)*sin(tout_c);
%Inner basic grid
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rp_c*cos(tin)+Rp_c*1i*sin(tin); 
tin_f = linspace(0,2*pi,N_f+1);
tin_f = tin_f(1:end-1)'; 
rbase_in_f = Rp_f*cos(tin_f)+Rp_f*1i*sin(tin_f); 

%Construct image grid
%will return only the basic outer grid, else refined outer grid 

rout = []; 
for k = 1:P
    rout = [rout; rbase_out_c+q(k)];
end
 
[~, ~, ~, rimage_vec, refine,pairs] = getEnhancedGrid(q, opt);

%Get pair basis
[Upf,Ypf,~,~,~,~] = getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt);


%Get one-body pseduo inverse blocks -- enough to do this for single body (everybody has the same coarse grid).
[UU,YY] = getSelfPseudo(1,rbase_in_c,rbase_out_c);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.rimage_vec = rimage_vec;
geom.opt = opt;
geom.rvec_out = rout;
geom.q = q;
geom.pairs = pairs;

basis = struct();
basis.U = UU;
basis.Y = YY;
basis.Upf = Upf;
basis.Ypf = Ypf;


%Visualise 1-body and pair-basis

%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,[],[],N_c, N_f,a_c,a_f,rads)

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
% Create new grid points, for which the accuracy of the solution is
% to be evaluated. 
rcheck_b = [];
n_bound = 2000; %check in a large number of points
t = linspace(0,2*pi,n_bound)';
for k = 1:P
    rcheck_b = [rcheck_b; q(k)+rads(k)*(cos(t)+1i*sin(t))];
end

%% SOLVE SYSTEM
% Build the matrix to inspect conditioning/eigenvalues if requested.
if debug
    x = zeros(2*length(rout),1);
    tic
    for k = 1:2*length(rout)
        k
        x(:) = 0; 
        x(k) = 1; 
        uu = matvec_2D_pairprecond_enhanced(x,geom,basis);
        CC(:,k) = uu;
    end
    toc
    figure(13);
    clf; 
    imagesc(log10(abs(CC)))
    colorbar
    title([solver_name ': log_{10} |CC|'],'interpreter','none')
    skeel(CC)

    figure(5)
    [V,D] = eig(CC);
    D = diag(D); 
    plot(real(D),imag(D),'ro')
    title([solver_name ': eigenvalues of CC'],'interpreter','none')

    [s,I] = mink(abs(D),3);
    Vsmall = V(:,I).*s';

    
    Mc = round(a_c*N_c);
    t = linspace(0,2*pi,Mc+1);
    t = t(1:end-1)'; 
    
    %% Visualise eigvecs of system matrix
    eigvec_vis = 0; 
    if eigvec_vis
        for i = 1:3
            for k = 1:P  
                Vpx = Vsmall((k-1)*Mc+1:k*Mc,i);
                Vpy = Vsmall((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
                Vp = abs(Vpx+1i*Vpy);
                figure(31)
                subplot(1,3,i)
                scatter3(real(q(k))+sin(t),imag(q(k))+cos(t),Vp,40,Vp,'filled');
                hold on
                colorbar
                view(0,90)
    
                figure(32)
                subplot(1,3,i)
                plot(t,Vp)
                hold on
    
    
                figure(27)
                subplot(1,3,i)
                plot(t,Vpx)
                hold on
    
                figure(28)
                subplot(1,3,i)
                plot(t,Vpy)
                hold on
    
                figure(29)
                subplot(1,3,i)
                scatter3(real(q(k))+sin(t),imag(q(k))+cos(t),Vpx,40,Vpx,'filled');
                hold on
                colorbar
                view(0,90)
    
                figure(30)
                subplot(1,3,i)
                scatter3(real(q(k))+sin(t),imag(q(k))+cos(t),Vpy,40,Vpy,'filled');
                hold on
                colorbar
                view(0,90)
    
    
    
            end

        end
    end
end
%% Experiment with left preconditioner based on deflation

if lr
    warning('Long range preconditioning not tested')
    rin_c = []; 
    for k = 1:P
        rin_c = [rin_c; rbase_in_c+q(k)];
    end
    [Sinv,Nx,Ny,Mx,Z,Y,db] = get_long_range_precond(q,rin_c,rout,opt);
    %tau_coarse1 = AN*Rinv*(AM'*fout); same thing
    opt.db = db;
    tau_coarse = getCoarseSource(fout,Sinv,Nx,Ny,Mx,Z,Y,db,P,opt.N_c,opt.a_c);
end


[tau,it,resvec,real_res] = helsing_gmres(@(x) matvec_2D_pairprecond_enhanced(x,geom,basis),fout,2*size(rout,1),maxit,gmres_tol,1,rout);

plot_gmres = true; 

%With Krylov precond, do something like
%[tau, e2, precond] = precond_gmres(@(x) matvec_2D_pairprecond_images(x,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rvec_out,q,UU,YY,pairs,Upf,Ypf,s), fout, zeros(2*size(rvec_out,1),1), 2*size(rvec_out,1), gmres_tol, precond,debug);
%it = length(e2); 

if plot_gmres
      figure()
%     semilogy(e2);
      semilogy(resvec);
      title('Convergence resistance with pair corr','interpreter','latex')

      %what's the resiudal?
    %  u = matvec_2D_pairprecond_images(tau,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,Upf,Ypf,s);
      u = matvec_2D_pairprecond_enhanced(tau,geom,basis);
%       figure()
%       semilogy(abs(u-fout))
%       title('Residual in solution')
end 


%% POSTPROCESS
[rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y] = getPairTransformationStokes(tau,geom,basis);

lambda =  [tau_stokes_x; tau_stokes_y];

%% Compute net forces and torque
% FT = [F_x; F_y; T_z] per particle, assembled from coarse, fine and
% pair-specific fine source strengths.
FT = zeros(3*P,1);

% Coarse source contribution.
Kc = getKmat2D(rbase_in_c,0);
for k = 1:P
    seg = (k-1)*3+1:3*k;
    idx = coarse_ind{k};
    FT(seg) = FT(seg) + Kc'*[tau_stokes_x(idx); tau_stokes_y(idx)];
end

% Fine and image-source contributions are stored particle-by-particle in
% the same order as in getPairTransformationStokes.
has_neigh = sort(unique(pairs(:)));
Kf = getKmat2D(rbase_in_f,0);
offset = P*N_c;

% Number of (extra) sources per particle for index bookkeeping.
n_im = zeros(P,1);
for row = 1:size(pairs,1)
    i1 = pairs(row,1);
    i2 = pairs(row,2);
    n_im(i1) = n_im(i1) + length(rimage_vec{i1,i2});
    n_im(i2) = n_im(i2) + length(rimage_vec{i2,i1});
end

for i = 1:length(has_neigh)
    k = has_neigh(i);
    seg = (k-1)*3+1:3*k;

    % Contribution from fine proxy ring (N_f points).
    fine_idx = offset + (1:N_f);
    offset = offset + N_f;
    FT(seg) = FT(seg) + Kf'*[tau_stokes_x(fine_idx); tau_stokes_y(fine_idx)];

    % Contribution from pair-specific extra/image points.
    nimk = n_im(k);
    if nimk > 0
        extra_idx = offset + (1:nimk);
        offset = offset + nimk;
        Kim = getKmat2D(rvec_in(extra_idx),q(k));
        FT(seg) = FT(seg) + Kim'*[tau_stokes_x(extra_idx); tau_stokes_y(extra_idx)];
    end
end

if offset ~= length(tau_stokes_x)
    error('Force postprocessing indexing mismatch in solve_res_precond_enhanced.');
end


%extract force and torque separately
for k = 1:P
    F(k) = FT((k-1)*3+1) + 1i*FT((k-1)*3+2);
    T(k) = FT(3*k); 
end


%% Do the evaluation of the flow in check points , FMM is applied.

ftest_b = getVelocityField(rvec_in,rcheck_b,...
    tau_stokes_x,tau_stokes_y);

%% Compute residual at boundary
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

%Some visualisation stuff... 
if visualise

    figure()
    plot(ftest_b)
    hold on
    plot([fb_x; fb_y])
    title('Velocity field rhs and lhs, res pair corr')

    %Visualise component by component
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
% 
%     sgtitle('Error on boundary','interpreter','latex')

    % Visualise the error for x and y together
    figure(10)
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
    title('Resistance pair corr')

    %% Visualise density
    figure()
    semilogy(abs(lambda))
    ylabel('Magnitude of MFS coefficients','interpreter','latex')
    title('Pair corr resistance')

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


function doPairBasisTest(Upf,Ypf,i,p2,q,U,Y,rbase_in_c,rbase_out_f,rpair_fine)
%Test pair basis with a smooth coarse density
warning('Check consistency with the matvec!')
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
        step1 = U{i}'*[tau_p1_x;tau_p1_y]; %here I assume x and y follow each other?
        tau_mapped = Y{i}*step1;
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
delta = 0.01; 
q = [0; 2+delta; (2+delta)*1i]; %center coordinates

%or, instead, three cireles in triangle - this DOES NOT converge to a good
%sol! 
% x = 1+delta/2;
% y = sqrt((2+delta)^2-(1+delta/2)^2);
% q = [0; 2+delta; x+1i*y];

U = [1 0; 0 0; 0 1]; %translational velocities
W = [1; 1; 1]; %angular velocities

%or, a pair of circles
% q = [0; 2+delta]; 
% U = [1 0; 0 0];
% W = [0; 0]; 

rads = [1; 1; 1]; 
visualise = 1; 
images = 1; 
delta_pair = 0.2; 
test = 1; % 1 or 2
%% compare to a solution with 1 body precond only
if test == 1
    %[FT,lambda,it1,gmres_tol,err1] = solve_2D_res(q,U,W,rads,images, visualise);
    rng(9);
    P = 3;
    delta = 0.05; %P = 5
    delta = 1e-2;
    %delta = 1;
    q = [0; 2+delta; 6];
    %q = grow_cluster(P,delta,2);
    side = 2 + delta;               % neighbor center distance
    R = side / (2*sin(pi/P));         % ring radius
    %q = R * exp(1i * (0:P-1).' * (2*pi/P));
   % q = [q; -6+1.5i; -2-4i]; P = P+2;
    %q = q([1,2,4],:); P = 3; 
    U = rand(P,2); W = rand(P,1); rads = ones(P,1);
    lr = 20; 
    lr = 0; 
    images = 1; 

    [FT1,lambda,it1,gmres_tol,err1] = solve_2D_res(q,U,W,rads,images, lr,visualise);
    [FT2,lambda,it2,gmres_tol,err2] = solve_res_precond_enhanced(q,U,W,rads,delta_pair,lr,visualise);
    [FT3,lambda,it3,gmres_tol,err3] = solve_res_precond_images(q,U,W,rads,delta_pair,lr,visualise);


    str = sprintf('Relative residual with 1-body precond: %1.2e vs 2-body: %1.2e\n Converging in %u resp % u iterations',err1,err2,it1,it2);
    disp(str)
    
    alignfigs;
else
    
    %% determine 2-way error (solve resistance followed by mob)
    lr = 0;
    [FT,lambda,it1,gmres_res, err_res] = solve_res_precond_images(q,U,W,rads,delta_pair,lr,visualise);
    F = [FT(1:3:end) FT(2:3:end)];
    T = FT(3:3:end); 
    [UW,lambdahat,it1,gmres_mob, err_mob] = solve_mob_precond_images(q,F,T,rads,delta_pair,visualise);
    Ures = [U W]';
    str = sprintf('Two way error is %1.3e, with resistance residual %1.3e and mobility residual %1.3e',norm(Ures(:)-UW),err_res,err_mob)
    alignfigs;
end
end
