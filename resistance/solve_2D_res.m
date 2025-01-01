function [F,T,err,iters,ftest,lambda,rcheck_dom] = solve_2D_res(q,u,omega,rads,image)
%% SOLVE_2D_RES(q,u,omega,rads,image) solves resistance problem with 1-body precond.
%Circles centered at q \in C^P with velocities u \in R^(P x 2), omega \in R
% R^P. Vector of radii rads. Image is a boolean to determine if image
% enhancement is to be used.

%% SET PARAMS
% GMRES PARAMS
maxit = 1600; 
gmres_tol = 1e-6; 
gmres_tol = 1e-10; %don't choose too large! It will impact the force resolution!

% Params to determine grid
a = 1.2; 
P = length(q); 
solve_xy = 1; %order system as x y x y x y if false
if image
    np = 60; %has been 60 all the time before
    %np = 200; 
    
    %np = 100; 
    tol = 1e-12;
    %tol = 1e-17;
    % np = 50; 
    %tol = 1e-8;
    wobbly = 0; %if other shape than circle, ignore for now.
else
    np = 200;
    np = 1000;
    np = 250; 
    
    %np = 200; 
    %np = 100; 
    %np = 60; 
    tol = 1e-14;
    %tol = 1e-12; 
    wobbly = 0; 
    %tol = 1e-10; % tol her[rout,rin,rimage_vec] = get2DImageGrid(Rg,image)e should probably be larger to allow for smaller coefficients, which in turn will improve on the GMRES convergence 
    %tol = 1e-14;
    %tol = 1e-5;
    %tol = 1e-9;
end

%% Set grid
%tol = 1e-12;
sep = (1/np)*log(1/tol);

%Rg = 0.9; /(a1+a2+d);
Rg = max([1-sep,0.01]); %radius of proxy surface
%Rg = 0.8;

%solve with Stresslets + Potential dipoles
s = [0 0 1 1];
mu = 1; %viscosity

visualise = 1; 

[rout,rin,rimage,nimage,pair_points] = get2DImageGrid(q,rads,Rg,a,np, image);

% rout = rout(1:end/2);
% rin = rin(1:end/2);
% rimage = rimage(1:end/2);
% nimage = nimage(1:end/2); 
% pair_points = pair_points(1,:); 
% q = q(1); 
% P = 1; 


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

% Create new grid points for checking residual, for which the accuracy of the solution is
% to be evaluated. 
rcheck_b = [];
n_bound = 803;

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
    np = 200; 
    sep = 0.2; %0.25 worked well
    for k = 1:size(q,1)
        s = wobblycurve(r0-a,a,w(mod(k,3)+1),nout,q(k));  
        theta = pi*rand(1); 
        %theta = pi/3;
        x = [cos(theta) -sin(theta); sin(theta) cos(theta)]*[real(s.x)'-real(q(k)); imag(s.x)'-imag(q(k))];
        rout = [rout; x(1,:)'+1i*x(2,:)'+q(k)];
        %create inner points
        s = wobblycurve(r0-a,a,w(mod(k,3)+1),np,q(k)); 
        
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

%% Construct boundary conditions
%Set boundary condition for MFS.
fout = [];
fb = []; 
start = 1;
for k = 1:P
   
   %fout = [fout; rhs_f(rout(nout*(k-1)+1:nout*k))];
   stop = sum(pair_points(1:k,2));
   if solve_xy
       rhs_f = @(x) [u(k,1)-omega(k)*(imag(x)-imag(q(k))); u(k,2)+omega(k)*(real(x-q(k)))]; 
       fout = [fout; rhs_f(rout(start:stop))];

       %test
       %K = getKmat2D(rout(start:stop),q(k));
       %uu = K*[u(k,:) omega(k)]';

   else
       rhs_f = @(x) [u(k,1)-omega(k)*(imag(x)-imag(q(k))) u(k,2)+omega(k)*(real(x-q(k)))]; 
       fout_temp = rhs_f(rout(start:stop))';
       fout = [fout; fout_temp(:)];

       fb_true = rhs_f(rcheck_b(n_bound*(k-1)+1:n_bound*k))';
       fb = [fb; fb_true(:)];
       
   end


   start = stop+1; 

end


%Construct chech boundaries
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



% Compute pseodoinverse for single particle

[Uii,Yii] = getSelfPseudo(rin,rout,rimage,nimage,P,pair_points,solve_xy,s);


%% Solve the system 

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


%x = rand(size(fout));

%res = matvec_2D_Stokes(x,rvec_in,rvec_out,rimage,nimage,q,UU,Y,pairs);

% Compute a solution with GMRES
%precond = KrylovPrecond();

%now, try to limit the number of iters for the preconditioner. Does
%that help?
% [x_gmres, e2, precond] = precond_gmres(@(x) matvec_2D_Stokes(x,rvec_in,rvec_out,rimage,nimage,q,UU,Y,pairs,s,wobbly), fout, zeros(2*size(rvec_out,1),1), 2*size(rvec_out,1), gmres_tol, precond,1);
%it = length(e2); 
%fprintf("iterations = %d\n", length(e2))
[x_gmres,iters,resvec,real_res] = helsing_gmres_mv(@(x) matvec_2D_Stokes(x,rin,rout,rimage,nimage,q,Uii,Yii,pair_points,s),fout,2*size(rout,1),maxit,gmres_tol,1,rout);
% [x,flag,relres,iter,resvec] = gmres(@(x) matvec_2D_Stokes(x,rvec_in,rvec_out,rimage,nimage,q,UU,Y,pairs,s,wobbly), fout,[],gmres_tol,500);


%delete precond

P = length(q); 
N_small = size(rin,1)/P;
PM = length(rout);

if solve_xy
    [tau_stokes_x,tau_stokes_y,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y] = getTransformedDensity(x_gmres,rimage,Uii,Yii,P,N_small,PM,pair_points,s);
    
    tau_proxy = [tau_stokes_x; tau_stokes_y];
    tau_image = [tau_stress_x; tau_stress_y; tau_pot_x; tau_pot_y];
    
    lambda = [tau_proxy; tau_image];
end

%And evaluate in new points rcheck_dom and rcheck_b

%% Do the evaluation of the flow in check points with fmm 
ftest_b = getFMMVelocity(rin,rcheck_b,tau_stokes_x,tau_stokes_y,rimage,nimage,tau_pot_x, tau_pot_y,tau_stress_x,tau_stress_y);
ftest = getFMMVelocity(rin,rcheck_dom,tau_stokes_x,tau_stokes_y,rimage,nimage,tau_pot_x, tau_pot_y,tau_stress_x,tau_stress_y);

%% Compute relative residual
fbound_x = ftest_b(1:length(rcheck_b));
fbound_y = ftest_b(length(rcheck_b)+1:end);


fb_x = [];
fb_y = [];
 

for k = 1:P
    rhs_f = @(x) [u(k,1)-omega(k)*(imag(x)-imag(q(k))); u(k,2)+omega(k)*(real(x-q(k)))]; 
    fb_true = rhs_f(rcheck_b(n_bound*(k-1)+1:n_bound*k));
    fb_x = [fb_x; fb_true(1:n_bound)];
    fb_y = [fb_y; fb_true(n_bound+1:end)];   
end
err = max(sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./max(sqrt(fb_x.^2+fb_y.^2))


%% Get force and torque 
for k = 1:P
    lambda_x = lambda(1:end/2);
    lambda_y = lambda(end/2+1:end);
    F1(k) = sum(tau_stokes_x((k-1)*np+1:k*np));
    F2(k) = sum(tau_stokes_y((k-1)*np+1:k*np));
    T(k) = sum(-tau_stokes_x((k-1)*np+1:k*np).*imag(rin((k-1)*np+1:k*np)-q(k)) + ...
      tau_stokes_y((k-1)*np+1:k*np).*real(rin((k-1)*np+1:k*np)-q(k)));
end
F = F1+1i*F2;




visualise = 0; 
if visualise
    figure(9)
    subplot(2,2,1)
    scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fb_x+fbound_x)),30,log10(abs(fb_x+fbound_x)),'filled')
    colorbar
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('error in x velocity')


    subplot(2,2,2)
    scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fb_y+fbound_y)),30,log10(abs(fb_y+fbound_y)),'filled')
    colorbar
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('error in y velocity')

    % Visualise the actual velocity
    subplot(2,2,3)
    scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fb_x)),30,log10(abs(fb_x)),'filled')
    colorbar
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('x velocity MFS')

    subplot(2,2,4)
    scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fbound_x)),30,log10(abs(fbound_x)),'filled')
    colorbar
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('x velocity reference')


    sgtitle('Error on boundary','interpreter','latex')

    % Visualise the error for x and y together
    figure(10)
    m = max(sqrt(fb_x.^2+fb_y.^2));
    scatter3(real(rcheck_b),imag(rcheck_b),log10((sqrt((fb_x+fbound_x).^2+(fb_y+fbound_y).^2))./m),...
        10,log10((sqrt((fb_x+fbound_x).^2+(fb_y+fbound_y).^2))./m),'filled')
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

    %% Visualise density
    lambda_x = lambda(1:length(rin));
    lambda_y = lambda(length(rin)+1:2*length(rin));

    lambda_image = lambda(2*length(rin)+1:end);

    figure(12)
    p = length(rin)/length(q);

    clf; 
    lambda_tot = vecnorm([lambda_y lambda_x],2,2);
    for k = 1:length(q)
        scatter(real(rin((k-1)*p+1:k*p)),imag(rin((k-1)*p+1:k*p)),...
            20,log10(abs(lambda_tot((k-1)*p+1:k*p))),'filled');
        hold on
    end
    %Visualise extra singularities at image points
    if size(rimage,1)
        i
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
    
end

end