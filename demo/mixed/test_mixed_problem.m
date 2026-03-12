%test combined mobility / resistance problem with one or multiple circles moving
%and two or four being fixed obstacles (the first Pw ojects). A unit force in the
%x direction is prescribed on the moving bodies. 

clear;
close all;

a = 4; %sets start and end position for the moving bodies.
delta_vec = logspace(-3,0,20); %set tightest separation between moving body and wall.
%delta_vec = 0.01;
P = 6; %number of objects (walls + particles)
images = 1; %solve accurately?
run_name = 'P6_accurate_obstacle_delta';
run_name = 'test';
%run_name = 'P6_crude_obstacle_delta';
%run_name = 'P6_more_crude_obstacle_delta';
%run_name = 'P6_accurate_obstacle_delta';
%run_name = 'P6_accurate_obstacle_delta';
dynamic_plot = 0; %draw every time step?
stream_lines = 0; %compute streeamlines?

%%
for k = 1:length(delta_vec)
    %delta = 0.01; 
    delta = delta_vec(k); 
    if k == 1
        dynamic_plot = 1; 
    else
        dynamic_plot = 0; 
    end

   
    %% Set start pos
    % Some different test configurations with 1-2 particles and 2-4 walls
    if P == 6
        y = [0 0 2*a 2*a -a -2*a; 
            2+delta -2-delta 2+delta -2-delta 0 0];
        Pw = 4;
        Pp = 2;
        F = [1 0; 1 0]; 
        T = [0; 0];
        ystat = [0 0 a a 0 -a; 2+delta -2-delta 2+delta -2-delta 0 0];
        ystat = [0 0 a a 0 -2-delta; 2+delta -2-delta 2+delta -2-delta 0 0];
        Pend = 5;
        max_x = 3*a; 
        
    elseif P == 4
        y = [0 0 -a -1.7*a; 2+delta -2-delta 0 0];
        Pw = 2;
        Pp = 2;
        F = [1 0; 1 0]; 
        T = [0; 0];
        ystat = [0 0 0 -4; 2+delta -2-delta 0 0];
        Pend = 3;
        max_x = a; 
    else
        y = [0 0 -a; 2+delta -2-delta 0];
        Pw = 2;
        Pp = 1; 
        F = [1 0]; 
        T = 0; 
        ystat = [0 0 0; 2+delta -2-delta 0];
        Pend = 3; 
        max_x = a; 
    end
    y = y(:);

     
    visualise = 1; %plot the configuration inside solve script?
   
    %Check accuracy level for a single static configuration (compare coarse and fine grid solution)  
    %solve_mixed solves a combined mobility and resistance problem for
    %rigid body velocities U on the moving particles and Fres forces and
    %torques on all walls. Here, rin, rimage, nimage and s are needed to
    %evaluate the flow field in post processing 
    getres = 1;
    [U,Fres,lambda_x,lambda_y,lambda_image,rin,rimage,nimage,relres] = solve_mixed(Pw,Pp,ystat(1:2:end)+1i*ystat(2:2:end),F,T,images,getres,visualise); 
    Uref = U(1); 

   
    %Determine dynamics -- solve in each time step 
    f = @(t,y) getVelocity(y,images,a,visualise,P,Pend,max_x);
    options =odeset('RelTol',1e-10);
    [t,y]=ode45(f,[0:0.5:1000],y,options); 
    
    %%
    %Only store data for which the dynamics is interesting (for small
    %enough x-coordinates on the particles, meaning they have not yet
    %passed the obstacles.
    ind = find(y(:,2*(Pend-1)+1)>max_x);
    if ~isempty(ind)
        t = t(1:ind(1));
        y = y(1:ind(1),:);
    end

    %% Determine flow field at each time stamp
    
    if dynamic_plot
    
        N = 100;
        warning('Cover everything with colloc points!')
        sx = linspace(-a-1,max_x,N);
        sy = linspace(-a,a,N); 
    
        figure(2)
    
    
        for i=1:length(t)
            [xx,yy] = meshgrid(sx,sy);
            coord = [xx(:)+1i*yy(:)];
            %coord = [linspace(-a,a,N)'+a*1i; linspace(-a,a,N)'-a*1i; 1i*linspace(-a,a,N)'+a; -1i*linspace(-a,a,N)'-a-1];
            q = y(i,1:2:end)+1i*y(i,2:2:end);
    
            indall = 1:length(coord);
            for j = 1:P
                ind = find(abs(coord-q(j))<1+1e-3);
                indall = setdiff(indall,ind);
            end
            
            %scatter(real(coord(indall)),imag(coord(indall)))
            rout = coord(indall);
            getres = 0; 
            [U,Fres,lambda_x,lambda_y,lambda_image,rin,rimage,nimage,~] = solve_mixed(Pw,Pp,q,F,T,images,getres,visualise);
            u = getField(t,[real(rout);imag(rout)],rin,rimage,nimage,lambda_x,lambda_y,lambda_image);
            uhist{i} = u;
            coordhist{i} = [real(rout);imag(rout)];
    
            % solve for streamlines?
            if stream_lines
                g = @(t,p) getField(t,p,rin,rimage,nimage,s,lambda);
                [tvec,streams] = ode45(g,[0,400],[real(rout);imag(rout)],options);
                %also, get velocities at these points
                for j = 1:length(tvec)
                    u = getField(tvec(j),[streams(j,1:end/2)';streams(j,end/2+1:end)'],rin,rimage,nimage,lambda_x,lambda_y,lambda_image);
                    uhist(:,j) = u;         
                end
            end
    
           
    
        end
    end

    % Save data to make comparisons later and also animate solution
    %save_str = sprintf('interp_accurate_obstacle_delta%1.2e.mat',delta);
    save_str = sprintf('%s%1.2e.mat',run_name,delta);
    if dynamic_plot
        save(save_str,'t','y','images','delta','uhist','coordhist','Uref','relres','a','Pw','Pp')
    else
        save(save_str,'t','y','images','delta','Uref','relres','a','Pw','Pp')
    end

end
%%
show_delta_dep = 1;
P = 6;

if show_delta_dep
    close all; 
    for k = 1:length(delta_vec)
        delta = delta_vec(k);
        save_str = sprintf('P%u_accurate_obstacle_delta%1.2e.mat',P,delta);
        load(save_str)
        tacc(k) = t(end);
        refval = Uref;
        res_acc(k) = relres;

        [maxmax,max_y,min_dist] = getMaxMax(y,t,Pw);
        maxmax_y_acc(k) = maxmax; 
        min_acc(k) = min(min_dist);
    
        save_str = sprintf('P%u_crude_obstacle_delta%1.2e.mat',P,delta);
        load(save_str)
        tcrude(k) = t(end);
        res_crude(k) = relres;
        Uerr(k) = abs(Uref-refval)./refval;
        [maxmax,max_y,min_dist] = getMaxMax(y,t,Pw);
        maxmax_y_crude(k) = maxmax; 
        min1(k) = min(min_dist); 
    
        save_str = sprintf('P%u_more_crude_obstacle_delta%1.2e.mat',P,delta);
        load(save_str)
        tcrude2(k) = t(end);
        res_crude2(k) = relres;
        [maxmax,max_y,min_dist] = getMaxMax(y,t,Pw);
        maxmax_y_crude2(k) = maxmax; 
        Uerr2(k) = abs(Uref-refval)./refval;
        min2(k) = min(min_dist);
    end
    
    figure()
    semilogx(delta_vec,tacc,'b+--')
    hold on
    semilogx(delta_vec,tcrude,'ro-')
    semilogx(delta_vec,tcrude2,'m*-')
    xlabel('$\delta$','interpreter','latex')
    ylabel('Time units','Interpreter','latex')
    grid on
    
    figure()
    loglog(delta_vec,abs(tacc-tcrude)./tacc,'ro-')
    hold on
    loglog(delta_vec,abs(tacc-tcrude2)./tacc,'m*-')
    xlabel('$\delta$','interpreter','latex')
    axis tight
    ylabel('Relative difference in end time','Interpreter','latex')
    grid on

    figure()
    loglog(delta_vec,res_acc,'b+--')
    hold on
    loglog(delta_vec,res_crude,'ro-')
    loglog(delta_vec,res_crude2,'m*-')
    xlabel('$\delta$','interpreter','latex')
    ylabel('Rel res at contact point','Interpreter','latex')
    grid on

    figure()
    loglog(delta_vec,Uerr,'ro-')
    hold on
    loglog(delta_vec,Uerr2,'m*-')
    xlabel('$\delta$','interpreter','latex')
    ylabel('Velocity error at contact point','Interpreter','latex')
    grid on

    figure()
    loglog(delta_vec,maxmax_y_acc,'b+--')
    hold on
    loglog(delta_vec,maxmax_y_crude,'ro-')
    loglog(delta_vec,maxmax_y_crude2,'m*-')
    xlabel('$\delta$','interpreter','latex')
    ylabel('max y coord of particle','Interpreter','latex')
    grid on

    figure()
    semilogx(delta_vec,min_acc,'b+--')
    hold on
    semilogx(delta_vec,min1,'ro-')
    semilogx(delta_vec,min2,'m*-')
    xlabel('$\delta$','interpreter','latex')
    ylabel('Smallest distance between objects','Interpreter','latex')
    grid on
    ylim([-0.1,0.1])


end
alignfigs;

%The RBM/force torque has been verified against a solution determined with a completed double layer BIE forumation

%% Animate the solution 
figure(1);
load('P6_accurate_obstacle_delta1.00e-03.mat')
for i = 1:length(t)
    u = uhist{i};
    maxval(i) = max(abs(u(1:end/2)+1i*u(end/2+1:end)));
end
Tend = t(end);
maxval = max(maxval); 

load('P6_more_crude_obstacle_delta1.00e-03.mat')
tt = linspace(0,2*pi,200);
%v = VideoWriter('Movie_P3_accurate.avi');  % can also use .mp4
v = VideoWriter('Movie_P6_more_crude.avi');  % can also use .mp4
v.FrameRate = 30;                % Set frames per second
v.Quality = 100; 
open(v);

for i = 2:length(t)
    clf
    for k = 1:Pw
        x_c = y(i,2*(k-1)+1);
        y_c = y(i,2*k);
        X = x_c+cos(tt);
        Y = y_c+sin(tt);
        fill(X,Y,[0 0 0])
        xlim([-a,max_x])
        %ylim([-a,a])
        axis equal

        hold on
    end
    for p = 1:Pp
        x_c = y(i,2*Pw+(p-1)*2+1);
        y_c = y(i,2*Pw+p*2);
        X = x_c+cos(tt);
        Y = y_c+sin(tt);
        fill(X,Y,'r')
    end
    plot([-1.2*a,max_x],[0,0],'k--','LineWidth',2)
    u = uhist{i};
    coord = coordhist{i};
    
    surfir(coord(1:end/2),coord(end/2+1:end),abs(u(1:end/2)+1i*u(end/2+1:end))-maxval);
    %caxis([-1.1*maxval,0])
    caxis([-maxval,-maxval/3]);
    
    coord_x = coord(1:end/2);
    coord_y = coord(end/2+1:end);
    u_x = u(1:end/2);
    u_y = u(end/2+1:end);
    c = colorbar;
    c.Ticks = [-maxval,-maxval/3];
    c.TickLabels = {'0','$2/3 \max \| \vec u\|$'};
    c.TickLabelInterpreter = 'latex';
    cc = 3;
    %cc = 5; %better choice? 
    x_coa = coord_x(1:cc:end);
    y_coa = coord_y(1:cc:end); 
    ux_coa = u_x(1:cc:end); 
    uy_coa = u_y(1:cc:end); 
    quiver(x_coa,y_coa,ux_coa,uy_coa,1,'k') 
    axis off
    str = sprintf('t = %1.3f',t(i)/Tend);
    title(str)

    frame = getframe(gcf);   % Capture the current figure
    writeVideo(v, frame);

    pause(0.0001/5)
   % pause(1)
   

    
end

close(v);



function rhs = getVelocity(y,images,a,visualise,P,Pend,max_x)
%get RBM velocities to update particle coordinates
if y(2*(Pend-1)+1)>max_x
    rhs = zeros(size(y));
    
else

    if P == 3
        F = [1 0]; 
        T = 0; %set force and torque on the particles
        Pp = 1;
        Pw = 2;
    elseif P == 4
        F = [1 0; 1 0];
        T = [0; 0]; 
        Pp = 2;
        Pw = 2;
    else
        Pw = 4;
        Pp = 2;
        F = [1 0; 1 0];
        T = [0; 0]; 
    end

    getres = 0;
    [U,Fres,lambda_x,lambda_y,lambda_image,rin,rimage,nimage,relres] = solve_mixed(Pw,Pp,y(1:2:end)+1i*y(2:2:end),F,T,images,getres,visualise); 
    %Determine velocity on particle together with forces/torques on the two fixed obstacles
    
    rhs = zeros(2*Pw,1);
    for k = 1:Pp
        rhs = [rhs; U((k-1)*3+1:(k-1)*3+2)'];
    end
end

    

end

%%%%%%%%%%%

function u = getField(t,rout,rin,rimage,nimage,stok_x,stok_y,lambda_image)
%Computes velocity field in the exterior of the particles! Probably better
%to do this with direct summation.

    rout = rout(1:end/2)+1i*rout(end/2+1:end); 

    rot = [];
    stress_x = lambda_image(1:end/4);
    stress_y = lambda_image(end/4+1:end/2);
    pot_x = lambda_image(end/2+1:3*end/4);
    pot_y = lambda_image(3*end/4+1:end); 
    u = getVelocityField(rin, rout, stok_x, stok_y, rimage, nimage, ...
                           rot, stress_x, stress_y,pot_x, pot_y);

end

function [maxmax,max_y,min_dist] = getMaxMax(y,t,Pw)
% % Check shortest distances
    for i = 1:length(t)
        yt = y(i,:);
        q = yt(1:2:end)+1i*yt(2:2:end);
        ycoord = yt(2*(Pw+1):2:end);
        max_y(i) = max(abs(ycoord));
        d = abs(q-q.');
%             for j = 1:P
% 
%                 d(j,:) = vecnorm(y(j,:)-q,2,2);
%             
%                 %Smallest dist for each particle?
%                 m = mink(d(j,:),2)-2;
%                 mind(j) = m(2);
%             end
        d = d(:);
        min_dist(i) = min(d(d>0))-2; 
    end
    maxmax = max(max_y); 
end