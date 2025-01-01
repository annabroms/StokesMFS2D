function [UW,lambda,it,err] = solve_2D_mob(q,F,T,rads,image)

%% PARAMS
%GMRES
maxit = 1600; 
gmres_tol = 1e-6;
gmres_tol = 1e-16; % not enough if we want to resolve forces

ver = 1; %version of the algorithm ver = 0 corresponds to the version where also image points are projected

%Params to set Rg and grid resolution
a = 1.2; 
if image
    np = 60; %has been 60 all the time before
    %np = 200; 
    np = 100; 
    
   % np = 100; 
    tol = 1e-12;
   % tol = 1e-16; 
    %tol = 1e-17;
    % np = 50; 
    %tol = 1e-8;
else
    np = 200;
    np = 1000;
    np = 100; 
    
    %np = 200; 
    %np = 100; 
    %np = 60; 
    tol = 1e-14;
    %tol = 1e-12; 
    %tol = 1e-10; % tol here should probably be larger to allow for smaller coefficients, which in turn will improve on the GMRES convergence 
    %tol = 1e-14;
end
%% SET GRID

%tol = 1e-12;
sep = (1/np)*log(1/tol);

%Rg = 0.9; /(a1+a2+d);
Rg = max([1-sep,0.01]); %radius of proxy surface
%Rg = 0.8;

%solve with Stresslets +Doublets at image points
s = [0 0 1 1];

visualise = 1; 

%create grids
%[rout,rin,rimage,nimage,pair_points,pairs,rimage_vec,refine,rin_base] = get2DImageGrid(q,rads,Rg_c,a_c,N_c,image,Rg_f,a_f,N_f,basic,delta_pair)
[rout,rin,rimage,nimage,pair_points] = get2DImageGrid(q,rads,Rg,a,np,image);
%[rout,rin,rimage,nimage,pair_points] = get2DImageGrid(q,rads,Rg, np, image);

% %debug
% rout = rout(1:end/2);
% rin = rin(1:end/2);
% rimage = rimage(1:end/2);
% nimage = nimage(1:end/2); 
% pair_points = pair_points(1,:); 
% q = q(1); 

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
visualise = 0; 
rcheck_b = [];
n_bound = 803;

% Grid to check residual
for k = 1:length(q)
    t = linspace(0,2*pi,n_bound)';
    rcheck_b = [rcheck_b; q(k)+rads(k)*(cos(t)+1i*sin(t))];
end

%% PREPARE PRECONDITIONING AND RHS

%get evaluation of lambda0
u = getRecompletionFlow(rin,rout,q,F,T); 

if ~isempty(rimage)
    [UU,Y,L] = getSelfPseudoMobility(rin,rout,q,rimage,nimage,length(q),pair_points,s,ver);
else
    [UU,Y,L] = getSelfPseudoMobility(rin(1:np)-q(1),rout(1:a*np)-q(1),0,rimage,nimage,1,pair_points,s,ver);
end


debug = 0;
if debug
    x = zeros(2*length(rout),1);
    tic
    for k = 1:2*length(rout)
        k
        x(:) = 0; 
        x(k) = 1; 
        uu = matvec_2D_mobility(x,rin,rout,rout,rimage,nimage,q,UU,Y,L,pair_points,s,1,ver);
        CC(:,k) = uu;
    end
    toc
    figure(14);
    clf; 
    imagesc(log10(abs(CC)))
    colorbar
    skeel(CC)
end


%% SOLVE
[x_gmres,it,resvec,real_res] = helsing_gmres_mv(@(x) matvec_2D_mobility(x,rin,rout,rout,rimage,nimage,q,UU,Y,L,pair_points,s,1,ver),u,2*size(rout,1),maxit,gmres_tol,1,rout);

figure()
semilogy(resvec)
title('GMRES convergence mobility, 1-body precond', 'Interpreter','latex')

P = size(q,1);
N_small = size(rin,1)/P;
PM = length(rout);

[tau_stokes_x,tau_stokes_y,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y] = getTransformedDensity(x_gmres,rimage,UU,Y,P,N_small,PM,pair_points,s,ver);

tau_proxy = [tau_stokes_x; tau_stokes_y];

tau_image = [tau_stress_x; tau_stress_y; tau_pot_x; tau_pot_y];
lambda = [tau_proxy; tau_image];

%get rigid body motion. 
UW = zeros(3*P,1); 
if ver
    K = getKmat2D(rin(1:N_small,:),q(1));
else
    start_ind = 0; 
end

for k = 1:P
    if ver      
        UW(3*(k-1)+1:3*k) = -K'*[tau_stokes_x((k-1)*N_small+1:k*N_small); tau_stokes_y((k-1)*N_small+1:k*N_small)];
    else
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
u_rhs = matvec_2D_mobility(x_gmres,rin,rout,rcheck_b,rimage,nimage,q,UU,Y,L,pair_points,s,0,ver);
S_0 = getRecompletionFlow(rin,rcheck_b,q,F,T); 
u_rhs = u_rhs-S_0; 

figure()
subplot(1,2,1)
plot(u_rhs);
hold on
plot(u_lhs);
subplot(1,2,2)
semilogy(abs(u_rhs-u_lhs'));


disp('Surface residual')
diff_vec = u_rhs-u_lhs';
max_abs = max(abs(u_rhs(1:end/2)+1i*u_rhs(end/2+1:end)));
err = max(abs(diff_vec(1:end/2)+1i*diff_vec(end/2+1:end)))/max_abs


visualise = 0; 
if visualise
    
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