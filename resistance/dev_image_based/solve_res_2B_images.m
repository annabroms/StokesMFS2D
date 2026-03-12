function [FT,lambda,it,gmres_tol,maxres] = solve_res_2B_images(q,U,W,rads,delta_pair,lr,visualise,gmres_tol,debug,gmres_verbose)
%SOLVE_RES_PRECOND_IMAGES Solves a 2D Stokes resistance problem with circular
%particles using MFS with 2-body preconditioning. To resolve
% challenging close interactions, a fine 2-body BVP is solved for fine
% source strengths on each close pair of particles. The fine sources,
% including sources at approximate image points, correct the representation
% obtained from a coarse grid, effectively preconditioning the system.
%
% Syntax:
%   [FT, lambda, it, gmres_tol, maxres] = solve_res_2B_images(q, U, W, rads, delta_pair, visualise)
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
%   along with image systems to improve solution accuracy for near-contact particle configurations.
%   The FMM is used for Stokeslet evaluation, while other source types use direct summation.
%
% Notes:
%   - This is an MFS generalisation of the Cheng-Greengard idea (1998/2000)
%   - Intended to match solve_mob_2B_images, but applied to a resistance problem rather than mobility.
%
% See also:
%   solve_res_1B              - 1-body preconditioned resistance solver
%   solve_res_peanut_images - 2-body preconditioner with peanut
%                              compression (in this version of the algorithm, only coarse sources
%                              are needed at the solve stage).
%   solve_mob_2B_images   - 2-body preconditioned mobility solver
%
% To test: Call without arguments.
%
% Anna Broms, April 9, 2025

if nargin==0, test_solve_res; 
    return; end

if nargin < 8 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 9 || isempty(debug), debug = false; end
if nargin < 10 || isempty(gmres_verbose), gmres_verbose = 0; end

P = length(q);

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
P = length(q); 

opt = get2Dparams();
opt.gmres_verbose = gmres_verbose;


%Play with N_c, N_f, a (a_f). 
N_c = 60;  %100 better here? 
%N_c = 80; 
%N_c = 150; 
%N_c = 100; 
N_f = 150; 
%N_f = 100;
%N_f = N_c; %debug

%N_c = 250; 
a_c = 1.2;
 
%a = 2; 
a_f = 1.2; %upsampling factor for the fine grid 

% Determine separation to proxy surface.
tol_c = 1e-12; %I think this works reasonably

%tol_c = 1e-10; %Curve moves closer to the surface -> smaller coeff 
%tol_c = 1e-16; %Curve moves further from surface -> larger coeff. 


s = [0 0 1 1 0 0 0]; %set type of singularities at image points
%s = [1 0 1 1]; %Other singularities? Currently not supported! But code can
%be changed!

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
 
opt.s = s;
opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c;
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 
opt.rads = rads;
opt.s = s; 
opt.P = P; 
opt.N_peanut = 0; 
opt.precomp = 1; %faster if evaluation of one body basis on fine grid is compted only once. 
% %Less storage required.
opt.pc = 1; %prepare grid to do pair corrections
opt.delta_pair = delta_pair; 
opt.lr = lr; %using long range preconditioning? 
opt.n_clusters = 30; 


%% CREATE GRID
%Outer basic grid
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';
rbase_out_c = rads(1)*cos(tout_c)+1i*rads(1)*sin(tout_c);
%Inner basic grid
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rp_c*cos(tin)+Rp_c*1i*sin(tin); 

%Construct image grid
%will return only the basic outer grid, else refined outer grid 
 
[rout,rin,rimage,~,~,pairs,rimage_vec,refine,rbase_in_f] = get2DImageGrid(q,rads,opt);
% pairs = [2,3];
% figure()
% for k = 1:P
%     plot(real(rbase_out_c+q(k)),imag(rbase_out_c+q(k)),'k-')
%     hold on   
%     scatter(real(rbase_in_c+q(k)),imag(rbase_in_c+q(k)),20,'filled','b')  
% end
% %axis off
% axis equal
% ylim([-5.5,3])
% set(gcf,'Color','white')
% 
% figure()
% for k = 1:P
%     if k ==2
%         fill(real(rbase_out_c+q(k)),imag(rbase_out_c+q(k)),[229 229 229]/255, 'EdgeColor', 'none');
%     end
%     hold on
% 
% 
%     if find(k==pairs(:))
%         scatter(real(rbase_in_f+q(k)),imag(rbase_in_f+q(k)),20,'filled','r') 
%         plot(real(rbase_out_c+q(k)),imag(rbase_out_c+q(k)),'k-','LineWidth',2)
%     else
%         plot(real(rbase_out_c+q(k)),imag(rbase_out_c+q(k)),'k-','LineWidth',1)
%         scatter(real(rbase_in_c+q(k)),imag(rbase_in_c+q(k)),20,'filled','b','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5) 
%         t_cc = linspace(0,2*pi,50);
%         rbase_out_cc = cos(t_cc)+1i*sin(t_cc);
%         plot(real(rbase_out_cc+q(k)),imag(rbase_out_cc+q(k)),'--','Color',[0.9 0.9,0.9])
%     end
% end
% for i = 1:size(pairs,1)
%     j = pairs(i,1);
%     k = pairs(i,2);  
%     scatter(real(rimage_vec{j,k}),imag(rimage_vec{j,k}),20,'filled','r');
%     scatter(real(rimage_vec{k,j}),imag(rimage_vec{k,j}),20,'filled','r');
% end
% %
% axis equal
% %ylim([-5.5,1.])
% %axis off
% set(gcf,'Color','white')


%Get pair basis
[Upf,Ypf,~,~,~,~,nimage] = getPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,[],[],[]);


%Get one-body pseduo inverse blocks -- enough to do this for single body (everybody has the same coarse grid).
[UU,YY] = getSelfPseudo(1,rbase_in_c,rbase_out_c);


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
if nargin < 6
    rcheck_dom = 100+100i; %some point far away
end

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
    ncols = 2*length(rout);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        uu = matvec_res_2B_images(x,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,Upf,Ypf,s);
        CC(:,k) = uu;
    end
    toc
    figure(13);
    clf; 
    imagesc(log10(abs(CC)))
    colorbar
    title([solver_name ': log_{10} |CC|'],'interpreter','none')
    cc = skeel(CC);
    fprintf('Estimated condition number of system matrix: %1.3e \n',cc);
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
%% Experiment with left preconditioner based on deflation

if lr
    rin_c = []; 
    for k = 1:P
        rin_c = [rin_c; rbase_in_c+q(k)];
    end
    [Sinv,Nx,Ny,Mx,Z,Y,db] = get_long_range_precond(q,rin_c,rout,opt);
    %tau_coarse1 = AN*Rinv*(AM'*fout); same thing
    opt.db = db;
    tau_coarse = getCoarseSource(fout,Sinv,Nx,Ny,Mx,Z,Y,db,P,opt.N_c,opt.a_c);
end




% Check older versions. Any symmetry speedups implemented? 

%res = matvec_res_2B_images(x,rbase_in_c,refine,rimage_vec,nimage,opt,rvec_out,q,U,Y,pairs,Upf,Ypf)
%res = matvec_2D_pairprecond(x,rvec_in,rvec_out,q,UU,Y,B);
%[tau,flag,relres,it,resvec2] = gmres(@(x) matvec_2D_pairprecond3(x,rbase_in_c,rbase_in_f,rbase_out_f,rvec_out,q,UU,YY,B,pairs,A,Uf,Yf,Ncf,Upf,Ypf),fout,[],gmres_tol,maxit);
disp(' == Solving... == ');
if lr
   Pf = applyPmat(fout,rin_c,rout,Sinv,Nx,Ny,Mx,Z,Y,opt);    
   [tau,it,resvec,real_res] = helsing_gmres(@(x) lr_matvec_2D_pairprecond(x,rin_c,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,Upf,Ypf,s,Sinv,Z,Y),Pf,2*size(rout,1),maxit,gmres_tol,opt,rout);   
else                                                                                 
    [tau,it,resvec,real_res] = helsing_gmres(@(x) matvec_res_2B_images(x,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,Upf,Ypf,s),fout,2*size(rout,1),maxit,gmres_tol,opt,rout);
end
plot_gmres = true; 

%With Krylov precond, do something like
%[tau, e2, precond] = precond_gmres(@(x) matvec_res_2B_images(x,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rvec_out,q,UU,YY,pairs,Upf,Ypf,s), fout, zeros(2*size(rvec_out,1),1), 2*size(rvec_out,1), gmres_tol, precond,debug);
%it = length(e2); 

if plot_gmres
      figure()
%     semilogy(e2);
      semilogy(resvec);
      title('Convergence resistance with pair corr','interpreter','latex')

      %what's the resiudal?
      u = matvec_res_2B_images(tau,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,Upf,Ypf,s);
%       figure()
%       semilogy(abs(u-fout))
%       title('Residual in solution')
end 


disp(' == Postprocessing == ');
%% POSTPROCESS
[rvec_in,rimage_in,nimage_in,coarse_ind,tau_stokes_x,tau_stokes_y,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y] = getPairTransformation(tau,rbase_in_c,rbase_in_f,refine,...
    rimage_vec,nimage,opt,rout,q,UU,YY,pairs,Upf,Ypf);

% if lr
%warning('Think about postprocessing here')
%      %tau_stokes = Qmat*[tau_stokes_x; tau_stokes_y];
%      tau_stokes = applyQmat([tau_stokes_x; tau_stokes_y],rin,rout,Rinv,AN,AM);
% 
%      tau_stokes_x = tau_stokes(1:end/2)+tau_coarse(1:end/2);
%      tau_stokes_y = tau_stokes(end/2+1:end)+tau_coarse(end/2+1:end);
%  end

rot = [];
lambda_image = [tau_stress_x; tau_stress_y; tau_pot_x; tau_pot_y; rot];
lambda_proxy = [tau_stokes_x; tau_stokes_y]; %This is the fine density
%And evaluate in new points rcheck_dom and rcheck_b
lambda = [lambda_proxy; lambda_image]; 


%% Compute forces
K = getKmat2D(rbase_in_c,0);
FT = zeros(3*P,1); 
for k= 1:P
    FT((k-1)*3+1:3*k) = K'*[tau_stokes_x((k-1)*N_c+1:k*N_c); tau_stokes_y((k-1)*N_c+1:k*N_c)];
end

K = getKmat2D(rbase_in_f,0);
has_neigh = sort(unique(pairs(:)));
for i = 1:length(has_neigh)
    k = has_neigh(i); 
    FT((k-1)*3+1:3*k) = FT((k-1)*3+1:3*k)+ K'*[tau_stokes_x((k-1)*N_f+1+P*N_c:k*N_f+P*N_c); 
        tau_stokes_y((k-1)*N_f+1+P*N_c:k*N_f+P*N_c)];
end


%% Do the evaluation of the flow in check points , FMM is applied.
ftest_b = getVelocityField(rvec_in,rcheck_b,...
    tau_stokes_x,tau_stokes_y,rimage_in,nimage_in,rot,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y);
%ftest = getVelocityField(rvec_in,rcheck_dom,...
%    tau_stokes_x,tau_stokes_y,rimage_in,nimage_in,rot,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y);


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




maxres = max(sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./max(sqrt(fb_x.^2+fb_y.^2));
fprintf('Max surf rel res at new nodes %.3e\n', maxres);

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
    semilogy(abs([lambda_proxy;lambda_image]))
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
delta = 0.001; 
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
    %[FT,lambda,it1,gmres_tol,err1] = solve_res_1B(q,U,W,rads,images, visualise);
    rng(9);
    P = 40;
    delta = 0.01; %P = 5
    %delta = 1;
    q = grow_cluster(P,delta,2);
    % side = 2 + delta;               % neighbor center distance
    % R = side / (2*sin(pi/P));         % ring radius
    % q = R * exp(1i * (0:P-1).' * (2*pi/P));
   % q = [q; -6+1.5i; -2-4i]; P = P+2;
    %q = q([1,2,4],:); P = 3; 
    U = rand(P,2); W = rand(P,1); rads = ones(P,1);
    lr = 20; 
    lr = 0; 
    images = 0; 
    [FT,lambda,it1,gmres_tol,err1] = solve_res_1B(q,U,W,rads,images, lr,visualise);
    [FT,lambda,it2,gmres_tol,err2] = solve_res_2B_images(q,U,W,rads,delta_pair,lr,visualise);
    
    str = sprintf('Relative residual with 1-body precond: %1.2e vs 2-body: %1.2e\n Converging in %u resp % u iterations',err1,err2,it1,it2);
    disp(str)
    
    alignfigs;
else
    
    %% determine 2-way error (solve resistance followed by mob)
    lr = 0;
    [FT,lambda,it1,gmres_res, err_res] = solve_res_2B_images(q,U,W,rads,delta_pair,lr,visualise);
    F = [FT(1:3:end) FT(2:3:end)];
    T = FT(3:3:end); 
    [UW,lambdahat,it1,gmres_mob, err_mob] = solve_mob_2B_images(q,F,T,rads,delta_pair,visualise);
    Ures = [U W]';
    str = sprintf('Two way error is %1.3e, with resistance residual %1.3e and mobility residual %1.3e',norm(Ures(:)-UW),err_res,err_mob)
    alignfigs;
end
end
