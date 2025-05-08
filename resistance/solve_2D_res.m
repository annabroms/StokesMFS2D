function [FT,lambda, it, gmres_tol, maxres] = solve_2D_res(q,U,W,rads,image,visualise)
%SOLVE_2D_RES Solves a 2D Stokes resistance problem with circular particles
% using 1-body preconditioned MFS.
%
% Syntax:
%   [FT, lambda, it, gmres_tol, maxres] = solve_2D_res(q, U, W, rads, image, visualise)
%
% Inputs:
%   q         - Vector of length P, complex-valued center coordinates for the particles
%   U         - Px2 matrix of translational velocities (columns: x and y components)
%   W         - Px1 column vector of angular velocities
%   rads      - Px1 vector of particle radii
%   image     - Logical flag (true/false): use image system for higher accuracy
%   visualise - Logical flag: plot the configuration and solution details
%
% Outputs:
%   FT         - 3P×1 vector of computed net forces and torques 
%   lambda  - Solution vector of source strengths
%   it         - Number of GMRES iterations required using GMRES with low stagnation control (J. Helsing).
%   gmres_tol  - Set GMRES tolerance
%   maxres     - Maximum relative residual in a test (non-collocation) set of boundary nodes
%
% See also:
%   solve_precond_images - 2-body preconditioning enhanced with images for
%                          resistance
%   solve_precond_peanut - 2-body preconditioning with peanut compression
%                          for resistance
%   solve_2D_mob         - 1-body preconditioned mobility formulation
%
% To test, call without arguments
%
% Anna Broms, April 4, 2025

if nargin==0, test_solve_res; 
    return; end

P = length(q);

%% Checks

assert(size(W,1)==P,'Wrong size of angular velocity vector')
assert(size(U,1)==P,'Wrong size of trans vel vector')
assert(size(U,2)==2,'Wrong size of trans vel vector, should contain x y coordinates')

if sum(rads)~=P
    warning('Double check radii') %in the design of image and collocation points, it's assumed all radii are 1.
end


%% SET PARAMS
% GMRES PARAMS
maxit = 1600;
gmres_tol = 1e-6;
gmres_tol = 1e-10; %don't choose too large! It will impact the force resolution!

% Params to determine grid
opt = get2Dparams(); 
a = opt.a_c;
opt.image = image;
solve_xy = 1; %order system as x y x y x y if false
if image
    Nc = 60; %has been 60 all the time before
    %Nc = 200;

    %Nc = 100;
    tol = 1e-12;
    %tol = 1e-17;
    % Nc = 50;
    %tol = 1e-8;
    wobbly = 0; %if other shape than circle, ignore for now.
else
    Nc = 200;
    Nc = 1000;
    Nc = 250;
    Nc = 60;
    tol = 1e-14;
    tol = 1e-12;
    wobbly = 0; %other particle shapes...
    %tol = 1e-10; % tol her should probably be larger to allow for smaller coefficients, which in turn will improve on the GMRES convergence
end


%tol = 1e-12;
sep = (1/Nc)*log(1/tol);

Rp = max([1-sep,0.01]); %radius of proxy surface
%Rp = 0.8;
opt.Rp_c = Rp;
opt.N_c = Nc;
opt.pc = 0; %no pair correction

%solve with Stresslets + Potential dipoles at image points
s = [0 0 1 1]; % s = [S R T D]
opt.s = s; 
mu = 1; %viscosity

%% GET GRIDS AND VISUALISE
%get grids
[rout,wout,rin,rimage,nimage,pair_points] = get2DImageGrid(q,rads,opt);

%Check coeff magnitude and residual if we solve with a single SVD.

% Nim = getImageKernels2D(rimage,nimage,rout,1,s);
% N = singleLayer(rin,rout,1);
%
% Ntot = [N Nim];
% [Ytot,Btot]  = getPseudoFactors(Ntot,1e-8,1);


if visualise
    figure()
    %Visualise particles
    plot(real(rin),imag(rin),'b.');
    hold on
    plot(real(rout),imag(rout),'k.');
    plot(real(rimage),imag(rimage),'b.');
    quiver(real(rimage),imag(rimage),real(nimage),imag(nimage))
    axis equal
end

% Create new grid points, in which the accuracy of the solution is
% to be evaluated.
rcheck_b = [];
n_bound = 1300;

if wobbly
    %ONLY WITHOUT IMAGE POINTS (ignore for now)
    w = [2 3 4];
    %w = [3 3 3];
    a = 0.15;
    a = 0.1;
    r0 = 1;
    nout = 200*2;
    tout = linspace(0,2*pi,nout+1);
    rout = [];
    rin = [];
    Nc = 200;
    sep = 0.2; %0.25 worked well
    for k = 1:size(q,1)
        s = wobblycurve(r0-a,a,w(mod(k,3)+1),nout,q(k));
        theta = pi*rand(1);
        %theta = pi/3;
        x = [cos(theta) -sin(theta); sin(theta) cos(theta)]*[real(s.x)'-real(q(k)); imag(s.x)'-imag(q(k))];
        rout = [rout; x(1,:)'+1i*x(2,:)'+q(k)];
        %create inner points
        s = wobblycurve(r0-a,a,w(mod(k,3)+1),Nc,q(k));

        y = [cos(theta) -sin(theta); sin(theta) cos(theta)]*[real(s.x)'-real(q(k)); imag(s.x)'-imag(q(k))];
        nx = [cos(theta) -sin(theta); sin(theta) cos(theta)]*[real(s.nx)'; imag(s.nx)'];
        y = y(1,:)'+1i*y(2,:)'+q(k);
        nx = nx(1,:)'+1i*nx(2,:)';
        x = y-sep*s.sp.*nx;
        rin = [rin; x];
        n_part(k) = nout;
        pair_points(k,2) = nout;

        s = wobblycurve(r0-a,a,w(mod(k,3)+1),n_bound,q(k));
        x = [cos(theta) -sin(theta); sin(theta) cos(theta)]*[real(s.x)'-real(q(k)); imag(s.x)'-imag(q(k))];
        rcheck_b = [rcheck_b; x(1,:)'+1i*x(2,:)'+q(k)];
    end


    if visualise
        figure()
        plot(real(rout),imag(rout),'b.');
        hold on
        plot(real(rin),imag(rin),'r.');
        axis equal
    end
else
    for k = 1:P
        t = linspace(0,2*pi,n_bound)';
        rcheck_b = [rcheck_b; q(k)+rads(k)*(cos(t)+1i*sin(t))];
    end
end

%Construct check boundaries
rcheck_dom = 100+100i;

tt = linspace(0,2*pi)';
rcheck_dom = 8*cos(tt)+8*1i*sin(tt);

% sizePlane = 500;
% rcheck_dom = [];
% L = 5;
% while size(rcheck_dom,1)<sizePlane
%     p = L*rand(1,2);
%     if ~any(vecnorm(q-[p 0],2,2)<1)
%         rcheck_dom = [rcheck_dom; p(1)+1i*p(2)];
%     end
% end

%% PREPARE PRECONDITIONING AND RHS
%Set boundary condition for MFS.
fout = [];
fb = [];
start = 1;
for k = 1:P

   %fout = [fout; rhs_f(rout(nout*(k-1)+1:nout*k))];
   stop = sum(pair_points(1:k,2));
   if solve_xy
       rhs_f = @(x) [U(k,1)-W(k)*(imag(x)-imag(q(k))); U(k,2)+W(k)*(real(x-q(k)))];
       fout = [fout; rhs_f(rout(start:stop))];

       %test
       %K = getKmat2D(rout(start:stop),q(k));
       %uu = K*[U(k,:) W(k)]';

   else
       rhs_f = @(x) [U(k,1)-W(k)*(imag(x)-imag(q(k))) U(k,2)+W(k)*(real(x-q(k)))];
       fout_temp = rhs_f(rout(start:stop))';
       fout = [fout; fout_temp(:)];

       fb_true = rhs_f(rcheck_b(n_bound*(k-1)+1:n_bound*k))';
       fb = [fb; fb_true(:)];

   end


   start = stop+1;

end

%If order unknowns with first x, then y instead of x y x y etc
if solve_xy
    %Rearrange the rhs
    start = 1;
    foutx = []; fouty = [];
    for k = 1:P
       stop = start+2*pair_points(k,2)-1;
       foutx = [foutx; fout(start:start+(stop-1-start)/2)];
       fouty = [fouty; fout(start+(stop-1-start)/2+1:stop)];


       start = stop+1;

    end
    fout = [foutx; fouty];
else
    error('See script solveStokes_withImagesFast');
end


% Compute pseodoinverse for single particle

[Uii,Yii] = getSelfPseudo(P,rin,rout,rimage,nimage,pair_points,s);

%% Solve the system




% Compute a solution using GMRES with Krylov preconditioning
%precond = KrylovPrecond();
%now, try to limit the number of iters for the preconditioner. Does
%that help?
% [x_gmres, e2, precond] = precond_gmres(@(x) matvec_2D_Stokes(x,rvec_in,rvec_out,rimage,nimage,q,UU,Y,pairs,s,wobbly), fout, zeros(2*size(rvec_out,1),1), 2*size(rvec_out,1), gmres_tol, precond,1);
%it = length(e2);
%fprintf("iterations = %d\n", length(e2))
%delete precond

% GMRES with low stagnation control (J.Helsing). 
[x_gmres,it,resvec,real_res] = helsing_gmres(@(x) matvec_2D_Stokes(x,rin,rout,rimage,nimage,q,Uii,Yii,pair_points,s),fout,2*size(rout,1),maxit,gmres_tol,1,rout);
% [x,flag,relres,iter,resvec] = gmres(@(x) matvec_2D_Stokes(x,rvec_in,rvec_out,rimage,nimage,q,UU,Y,pairs,s,wobbly), fout,[],gmres_tol,500);

% Decay of residual with iteration number
figure()
semilogy(resvec)
title('GMRES convergence resistance, 1-body precond', 'Interpreter','latex')
xlabel('Iteration number','interpreter','latex')

PM = length(rout);

if solve_xy
    [tau_stokes_x,tau_stokes_y,rot,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y] = getTransformedDensity(x_gmres,rimage,Uii,Yii,P,Nc,PM,pair_points,s);

    tau_proxy = [tau_stokes_x; tau_stokes_y];
    tau_image = [rot;tau_stress_x; tau_stress_y; tau_pot_x; tau_pot_y];

    lambda = [tau_proxy; tau_image];
end


% Does i matter if we rescale the least squares problem? 
Nio = singleLayer(rin,rout,mu); 
Nimage = getImageKernels2D(rimage,nimage,rout,mu,s); 
   
Nii = [Nio Nimage];

w2 = vecnorm(Nii);
W2 = diag(1./w2);

W1 = blkdiag(diag(sqrt(wout)),diag(sqrt(wout)));

%[Y,Bi]  = getPseudoFactors(W1*Nii*W2,tol,visualise);
[Y,Bi]  = getPseudoFactors(Nii*W2,tol,visualise);
lambda2 = W2*Y*(Bi'*fout);

[Y,Bi]  = getPseudoFactors(Nii,tol,visualise);
lambda3 = Y*(Bi'*fout);

[Y,Bi]  = getPseudoFactors(W1*Nii, tol,visualise);
lambda4 = Y*(Bi'*W1*fout);

[Y,Bi]  = getPseudoFactors(W1*Nii*W2,tol,visualise);
lambda5 = W2*Y*(Bi'*W1*fout);

lambda = lambda4;
tau_stokes_x = lambda(1:opt.N_c*P);
tau_stokes_y = lambda(opt.N_c*P+1:2*opt.N_c*P);
tau_stress_x = lambda(2*opt.N_c*P+1:2*opt.N_c*P+size(rimage,1));
tau_stress_y = lambda(2*opt.N_c*P+size(rimage,1)+1:2*opt.N_c*P+2*size(rimage,1));
tau_pot_x =  lambda(2*opt.N_c*P+2*size(rimage,1)+1:2*opt.N_c*P+3*size(rimage,1));
tau_pot_y =  lambda(2*opt.N_c*P+3*size(rimage,1)+1:2*opt.N_c*P+4*size(rimage,1));



%And evaluate in new points rcheck_dom and rcheck_b

%% Do the evaluation of the flow in check points  rot,stress_x,stress_y,pot_x,pot_y
ftest_b = getVelocityField(rin,rcheck_b,tau_stokes_x,tau_stokes_y,rimage,nimage,rot,tau_stress_x,tau_stress_y,tau_pot_x, tau_pot_y);
ftest = getVelocityField(rin,rcheck_dom,tau_stokes_x,tau_stokes_y,rimage,nimage,rot,tau_stress_x,tau_stress_y,tau_pot_x, tau_pot_y);

%% Compute relative residual
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


%% Get force and torque
Kin = getKmat2D(rin(1:Nc),q(1));

for k = 1:P
    lambda_x = tau_stokes_x((k-1)*Nc+1:k*Nc);
    lambda_y = tau_stokes_y((k-1)*Nc+1:k*Nc);
%     F1(k) = sum(tau_stokes_x((k-1)*Nc+1:k*Nc));
%     F2(k) = sum(tau_stokes_y((k-1)*Nc+1:k*Nc));
%     T(k) = sum(-tau_stokes_x((k-1)*Nc+1:k*Nc).*imag(rin((k-1)*Nc+1:k*Nc)-q(k)) + ...
%         tau_stokes_y((k-1)*Nc+1:k*Nc).*real(rin((k-1)*Nc+1:k*Nc)-q(k)));
    
    %Same thing:
    FT((k-1)*3+1:k*3) = Kin'*[lambda_x; lambda_y];

end
%F = F1+1i*F2;

if visualise
    %visualise velocity and error component wise
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
%     sgtitle('Error on boundary','interpreter','latex')

    % Visualise the error for x and y together
    figure()
    m = max(sqrt(fb_x.^2+fb_y.^2));
    scatter3(real(rcheck_b),imag(rcheck_b),log10((sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./m),...
        10,log10((sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./m),'filled')
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    c = colorbar;
    ylabel(c,'Relative residual (log10)','Interpreter','latex')
    set(c,'TickLabelInterpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    set(gcf,'color','w');
    axis off
    title('1 body precond resistance')

    %% Visualise density
    lambda_x = lambda(1:length(rin));
    lambda_y = lambda(length(rin)+1:2*length(rin));

    lambda_image = lambda(2*length(rin)+1:end);

    figure()
    clf;
    lambda_tot = vecnorm([lambda_y lambda_x],2,2);
    for k = 1:length(q)
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
    title('1-body precond resistance')

end

end

function test_solve_res

close all; 
rng(6)
images = 0; %images not needed for well separated particles
q = [0; 2.01]; %center coordinates
U = [1 0; 0 0]; %translational velocities 
W = [1; 1]; %angular velocities 
rads = [1; 1]; 
visualise = 1; 
%[FT,lambda, it, gmres_tol, err] = solve_2D_res(q,U,W,rads,images,visualise);

%compare to a solution with image enhancement
images = 1; 
[FT,lambda, it, gmres_tol, err_im] = solve_2D_res(q,U,W,rads,images,visualise);

str = sprintf('Relative residual with image enhancement: %1.2e vs without: %1.2e',err_im,err);
disp(str)
alignfigs; 



end
