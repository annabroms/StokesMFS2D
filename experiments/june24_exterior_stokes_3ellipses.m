%% Setup a problem with three ellipses, similar to the june_17_exterior_3ellipse script for Laplace
close all;
clear;

svec = [1 0 0 1 1]; % set of singularities at clustered nodes [S R Tr D,T]: stokeslets, 
% rotlets, stresslet with random "normal" direction, potential dipoles and stresslets with two directions ex,ey

% Set particle configuration
delta = 1e-5; %not sure atm what the resulting particle-particle distance is
avec = [1 1 3]; %set axis
bvec = [2,2,1]; %set axis
taus = [0 0 -pi/25]; %set rotation
cvec = [0; 2+delta-0.2*1i; 1-3.07*1i];


r = 0.8; %sets proxy curve "radius", i.e. 1-r is the imaginary shift of the parameterized boundary curve
r = 0.7; 
N1 = 120; 
osf = 1.2;

Nclust = 30; %Set nbr of clustered nodes per particle per pair
Mclust = Nclust*5; % Possible to pick more clustered collocation points
threshold = 1e-6; 
threshold = 1e-1; %for clustered nodes: seems to work well! 


trunc = 1e-15; 
scalenorms = 1; 

%% set bc
P = 3; % number of particles
Up = [0 0; 1 0; 0 1]; Wp = [1; 0; 0]; % translational and angular velocities

M1 = ceil(osf*N1);
t = linspace(0, 2*pi, M1+1)';
t = t(1:end-1);


coll = [];
rcheck = [];
ntest = 2e4; 
tcheck = linspace(0,2*pi,ntest+1);
tcheck = tcheck(1:end-1)'; 

for i = 1:P
    zs = @(t) cvec(i)+(avec(i)*cos(t)+bvec(i)*1i*sin(t))*exp(taus(i)*1i);
    b_funcs{i} = zs;
    coll = [coll; zs(t)];
    rcheck = [rcheck; zs(tcheck)];

end

%% create rhs vector
start = 1; 

u_ref = zeros(P*ntest*2,1); 

fout = zeros(size(coll,1)*2,1); 
for k = 1:P
   stop = k*M1;
   rhs_f = @(x) [Up(k,1)-Wp(k)*(imag(x)-imag(cvec(k))); Up(k,2)+Wp(k)*(real(x-cvec(k)))];
   rhs_k = rhs_f(coll(start:stop));
   fout(start:stop) = rhs_k(1:end/2);
   fout(P*M1+start:P*M1+stop) = rhs_k(end/2+1:end);

   start = stop+1;
   rcheck_k = rcheck((k-1)*ntest+1:k*ntest);
   Ktest = getKmat2D(rcheck_k,cvec(k));
   u = Ktest*[Up(k,:)'; Wp(k)];
   u_ref((k-1)*ntest+1:k*ntest) = u(1:end/2);
   u_ref((k-1)*ntest+P*ntest+1:k*ntest+P*ntest) = u(end/2+1:end);  


end

%% Source points
%Use imaginary shift of the boundary curve to obtain source points, Liu style. 
t = linspace(0,2*pi,N1+1);
t = t(1:end-1)';
centers = [];
for i = 1:3
    % We need f for later assigning fine collocation points
    [~,~,~,f,L] = distributePoints(avec(i),bvec(i),taus(i),N1);
    flist{i} = f;
    Llist(i) = L;
    zs_i = b_funcs{i};
    centers_i = zs_i(t+(1-r)*1i);
    centers = [centers; centers_i];
end


%% visualise
figure()
plot(real(centers),imag(centers),'*')
hold on
%plot(real(centers1_test),imag(centers1_test),'ro')
plot(real(coll),imag(coll),'*')
axis equal
title('basic setup')

%% Solve the problem without clustering

[x_out,reserr_c,coefnorm_c,solres] = solve_dense_stokes(centers,coll,[],[],svec,fout,trunc,scalenorms);

u = eval_Stokes_solution(x_out,centers,rcheck,[],[],svec,N1,P);
diff_vec = u-u_ref;
diff_vec_c = diff_vec(1:end/2)+1i*diff_vec(end/2+1:end);
err_c = max(abs(diff_vec_c));


fprintf('Without clustering: Residual (collocation nodes) %1.3e, residual (new nodes) %1.3e and coefficient norm %1.3e\n', reserr_c,err_c,coefnorm_c)


%% Solution with accumulation points and clustering

pairs = [1 2; 2 3; 1 3]; % all neighbours
clusters = [];
coll_clust = coll; %start with basic grid
fout_x = fout(1:end/2);
fout_y = fout(end/2+1:end); 

%tclust = 0.5*logspace(log10(threshold),0,Nclust);
tclust = 0.8*logspace(log10(threshold),0,Nclust);

tcollclust = logspace(log10(1e-6),0,Mclust)';

for k = 1:size(pairs,1)
    inds = pairs(k,:);
    % get discretization for the two particles in the pair (we have assumed
    % equally many source and collocation points for everyone).
    z1 = coll((inds(1)-1)*M1+1:inds(1)*M1);
    z2 = coll((inds(2)-1)*M1+1:inds(2)*M1);

    % Get clustered nodes
    schwarz1 = aaa(conj(z1), z1);
    schwarz2 = aaa(conj(z2), z2);

    c1 = cvec(inds(1)); c2 = cvec(inds(2));

    [zacc1,zvec1,itr1] = compute_reflections(c2, schwarz1, schwarz2);
    [zacc2,zvec2,itr2] = compute_reflections(c1, schwarz2, schwarz1);

    cent_clust1 = zacc1-tclust'*(zacc1-c1)./abs(zacc1-c1)+tclust(1)*(zacc1-c1)./abs(zacc1-c1); %the clustering should end at the accumulation point! 
    cent_clust2 = zacc2-tclust'*(zacc2-c2)./abs(zacc2-c2)+tclust(1)*(zacc2-c2)./abs(zacc2-c2);

    %works equally well: 
    %cent_clust1 = zacc1+tclust'*(zacc1-zacc2)./abs(zacc1-zacc2)-tclust(1)*(zacc1-zacc2)./abs(zacc1-zacc2); %the clustering should end at the accumulation point! 
    %cent_clust2 = zacc2+tclust'*(zacc2-zacc1)./abs(zacc2-zacc1)-tclust(1)*(zacc2-zacc1)./abs(zacc2-zacc1);

    clusters = [clusters; cent_clust1; cent_clust2];

    %% Also get additional collocation points
    % To assign clustered nodes at the boundary, first find closest parameter value on the two curves
    % Objective function: squared distance to zacc1
    f1 = flist{inds(1)};
    f2 = flist{inds(2)};
    L1 = Llist(inds(1)); L2 = Llist(inds(2));

    dist1 = @(t) abs(f1(t)+c1 - zacc1).^2;
    
    % Minimize over t in [0, L]
    t_opt1 = fminbnd(dist1, -0.2, L1+0.2); %make sure to find minimum
    
    % Closest point on the curve
    z_closest1 = f1(t_opt1)+c1;
    tclust1 = mod([tcollclust; -tcollclust]+t_opt1,L1);
    
    %repeat for the other obstacle
    dist2 = @(t) abs(f2(t)+c2 - zacc2).^2;
    t_opt2 = fminbnd(dist2, -0.2, L2+0.2);
    z_closest2 = f2(t_opt2)+c2;
    tclust2 = mod([tcollclust; -tcollclust]+t_opt2,L2);
    
    coll_clust1 = c1 + f1(tclust1);
    coll_clust2 = c2 + f2(tclust2);

    coll_clust = [coll_clust; coll_clust1; coll_clust2];

    % Add rhs vector 
    i = inds(1);
    rhs1 = @(x) [Up(i,1)-Wp(i)*(imag(x)-imag(cvec(i))); Up(i,2)+Wp(i)*(real(x-cvec(i)))];
    rhs_1 = rhs1(coll_clust1);
    fout_x = [fout_x; rhs_1(1:end/2)];
    fout_y = [fout_y; rhs_1(end/2+1:end)];
    
    i = inds(2);
    rhs2 = @(x) [Up(i,1)-Wp(i)*(imag(x)-imag(cvec(i))); Up(i,2)+Wp(i)*(real(x-cvec(i)))];
    rhs_2 = rhs2(coll_clust2);
    fout_x = [fout_x; rhs_2(1:end/2)];
    fout_y = [fout_y; rhs_2(end/2+1:end)];

end

fout_fine = [fout_x; fout_y];


figure()
plot(real(centers),imag(centers),'*')
hold
plot(clusters, 'rx')
hold on
axis equal
plot(coll_clust, 'b.')
title('Enhanced setup')
hold off

%%

if svec(5)
    dir = [ones(length(clusters),1);1i*ones(length(clusters),1)];
else
    dir = rand(length(clusters),1);
end


[x_f,reserr_f,coefnorm_f,solres] = solve_dense_stokes(centers,coll_clust,clusters,dir,svec,fout_fine,trunc,scalenorms);

u = eval_Stokes_solution(x_f,centers,rcheck,clusters,dir,svec,N1,P);
diff_vec = u-u_ref;
diff_vec_f = diff_vec(1:end/2)+1i*diff_vec(end/2+1:end);
err_f = max(abs(diff_vec_f));


fprintf('Without clustering: Residual (collocation nodes) %1.3e, residual (new nodes) %1.3e and coefficient norm %1.3e\n', reserr_f,err_f,coefnorm_f)



%%
figure()
subplot(1,2,1)
scatter(real(rcheck),imag(rcheck),5,log10(abs(diff_vec_c)))
axis equal
colorbar
title('Basic setup')
sgtitle('Surface residual in new points')
clim([-9,0])
subplot(1,2,2)
scatter(real(rcheck),imag(rcheck),5,log10(abs(diff_vec_f)))
axis equal
colorbar
title('Enhanced setup')
sgtitle('Surface residual in new points')
clim([-9,0])

%% Visualise solution field
sizePlane = 20000;

%Generate random points in the exterior domain
L = 1.8;
rcheck_dom = [];
while size(rcheck_dom,1)<sizePlane
    p = -1.5*L+4*L*rand(1)+1i*(-2.5*L+4*L*rand(1));
    next = 0;
    for k = 1:P
        if next
            break;
        end
        inside = is_inside_rotated_ellipse(p, cvec(k), avec(k), bvec(k), taus(k));
        if inside
            next = 1;
        end
    end
    if ~next
        rcheck_dom = [rcheck_dom; p];
    end
end

figure()
subplot(1,2,1)
u = eval_Stokes_solution(x_f,centers,rcheck_dom,clusters,dir,svec,N1,P);
surfir(real(rcheck_dom),imag(rcheck_dom),u(1:end/2));
view(0,90)
colorbar
clim([-1,1])
axis square
axis tight
removePatches(cvec,0,avec,bvec,taus,1e6)  %only visualise the exterior field
title('$u_x$','Interpreter','latex')

subplot(1,2,2)
surfir(real(rcheck_dom),imag(rcheck_dom),u(end/2+1:end));
view(0,90)
colorbar
clim([-1,1])
axis square
axis tight
removePatches(cvec,0,avec,bvec,taus,1e6)  %only visualise the exterior field
title('$u_y$','Interpreter','latex')


sgtitle('Solution field')


alignfigs();

function removePatches(q,h,avec,bvec,taus,maxval) 
    t = linspace(0,2*pi,200);
    for k = 1:length(q)
        r_range = linspace(0,1+h,2);
        [R,T] = meshgrid(r_range,t);
        X = avec(k)*R.*cos(T);
        Y = bvec(k)*R.*sin(T);
        Z = (X+1i*Y)*exp(1i*taus(k));
        X = real(Z);
        Y = imag(Z); 

        patch(X(:)+real(q(k)),Y(:)+imag(q(k)),maxval*ones(size(X(:))),'w','EdgeColor', 'w');
    end
end

function inside = is_inside_rotated_ellipse(zp, c, a, b, tau)
    % zp  - complex point to test
    % c   - complex center of ellipse
    % a   - semi-major axis
    % b   - semi-minor axis
    % tau - rotation 
    
    % Translate and rotate zp into the local frame of the ellipse
    z_local = (zp - c) * exp(-1i * tau);
    
    % Extract real and imaginary parts as x' and y'
    x = real(z_local);
    y = imag(z_local);
    
    % Evaluate ellipse inequality
    inside = (x/a)^2 + (y/b)^2 <= 1;
end