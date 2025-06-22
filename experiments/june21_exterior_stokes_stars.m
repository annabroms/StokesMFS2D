%% Setup a problem with two stars
close all;
clear;

% Change svec for different singularities at the clustered nodes 

% the obstacles are spaced delta aside
delta = 1e-5;
%delta = 1; 
trunc = 1e-13; %SVD truncation level
scalenorms = 1; % do column scaling

Nclust = 100; %Set nbr of additional sources in each star clustered towards 
% image accumulation points

%Nclust = 200; 

fprintf('------Stokes resistance demo with two stars delta = %1.2e apart------\n',delta)
svec = [1 0 0 1 1]; % S R Tr D T at clustered nodes
%svec = [0 0 0 1 1];
P = 2; 

% Set rigid body motion
Up = [1 0; -1 0]; 
Wp = [1;  0]; 
%Wp = [0; 0];

%% Set geometry
%number of arms
l = 5;
%l = 7;

%amplitude of arms
a = 0.1;
%a = 0.3;

vars(1) = a;
vars(2) = l; 

rotate = 0; %different performance depending on geometry! Why this large difference?

c1 = 0;
if rotate

    % rotation
    tau1 =-pi/2/l; % with this choice, the arms will meet 
    tau2 = pi/2/l; 
    c2 = 2+2*a+delta; %stars delta apart

else

    tau1 = 0; 
    tau2 = 0;  
    f = @(t) (1 + a * sin(l * t)) .* cos(t);
    tmax = fminbnd(@(t) -f(t), 0, 2*pi);  % Find largest x val
    fmax = f(tmax);
    c2 = 2*fmax+delta; %Now, the particles are agian delta apart.
end

q = [c1; c2];


t = linspace(0, 2*pi, 401);
t = linspace(0, 2*pi, 801); %change here if a different number of arms
s1 = @(t) (1+a*sin(l*t)).*exp(1i*t);
z1 = c1+s1(t)*exp(1i*tau1);
z2 = c2+s1(t)*exp(1i*tau2); 

%% Make a set of source points
%r = 0.9;
%r = 0.75
N1 = 150;
%N1 = 100; 
N1 = 90; 

% Use the natural parameterization 
[zs1,n1,zp1,zpp1] = getPoints(a,l,tau1,N1,'star');
r = 0.7; % sets offset from boundary
%r = 0.8;
[~,~,~,f1,L] = distributePoints(a,l,tau1,N1,'star');
[~,~,~,f2,~] = distributePoints(a,l,tau2,N1,'star');


% Obtain source points via imaginary shift of boundary parameterization
t = linspace(0,2*pi,N1+1)';
t = t(1:end-1); 
fz = @(t) (1+a*sin(l*t)).*exp(1i*t)*exp(1i*tau1);
centers1 = fz(t+(1-r)*1i)+c1; 

fz = @(t) (1+a*sin(l*t)).*exp(1i*t)*exp(1i*tau2);
centers2 = fz(t+(1-r)*1i)+c2;
    

centers = [centers1; centers2];

%% Make a bunch of collocation points
osf = 1.2; %was 10
%osf = 10; 
M1 = N1 * osf; 

[z,n1,zp1,zpp1] = getPoints(a,l,tau1,M1,'star');
coll1 = z;
[z,n2,zp1,zpp1] = getPoints(a,l,tau2,M1,'star');
coll2 = z;


coll = [coll1+c1; coll2+c2];


%% Solution without accumulation points

M = length(coll);

% First, get rhs
bvec = [];
Mk = M/P;
start = 1; 
for k = 1:P
   stop = Mk*k;
   rhs_f = @(x) [Up(k,1)-Wp(k)*(imag(x)-imag(q(k))); Up(k,2)+Wp(k)*(real(x-q(k)))];
   rhs_k = rhs_f(coll(start:stop));
   bvec(start:stop) = rhs_k(1:end/2);
   bvec(M+start:M+stop) = rhs_k(end/2+1:end);  
   start = stop+1; 
end
bvec = bvec.';


%% solve with basic grid

[x_b,reserr_p,coefnorm_p] = solve_dense_stokes(centers,coll,[],[],[0 0 0 0 0],bvec,trunc,scalenorms);

%% Visualise geometry
figure()
plot(centers, 'rx')
hold on
%plot(centers1_test,'co')
axis equal
plot(coll, 'b.')
quiver(real(coll),imag(coll),real([n1;n2]),imag([n1;n2])) %visualise normals
title('Basic setup')
%% Visualise source strengths
figure()
scatter(real(centers),imag(centers),20,log10(abs(x_b(1:end/2)+1i*x_b(end/2+1:end))),'filled')
axis equal
colorbar
title('Source strength magnitude log10, basic')

%% Visualise boundary data
figure()
subplot(1,2,1)
scatter(real(coll),imag(coll),10,bvec(1:end/2),'filled')
colorbar
axis equal

subplot(1,2,2)
scatter(real(coll),imag(coll),10,bvec(end/2+1:end),'filled')
colorbar
axis equal


%% Solution with accumulation points and clustering

[schwarz1,pol1, res, zer, zj, fj, wj, errvec, wt] = aaa(conj(z1), z1);
[schwarz2,pol2, res, zer, zj, fj, wj, errvec, wt] = aaa(conj(z2), z2);
figure(1)
hold on
plot(real(pol1),imag(pol1),'k.')
plot(real(pol2),imag(pol2),'c.')



[zacc1,zvec1,itr1] = compute_reflections(c2, schwarz1, schwarz2);
[zacc2,zvec2,itr2] = compute_reflections(c1, schwarz2, schwarz1);
%[zacc1_2,~,itr3] = compute_reflections(c2+0.01, schwarz1, schwarz2); %converges to the same points
%[zacc2_2,~,itr4] = compute_reflections(c1+0.01, schwarz2, schwarz1);


threshold = 1e-10;

tclust = 0.5*logspace(log10(threshold),0,Nclust);
%tclust = 0.4*exp(-2*pi*(0:Nclust-1)/sqrt(Nclust)); %differently clustered nodes:  see "Resolution of sing. by rational functions"
%tclust = (zacc1-r)*0.8*exp(-2*sqrt(2)*pi*(sqrt(Nclust)-sqrt(1:Nclust))); %differently clustered nodes

cent_clust1 = [zacc1-tclust'*(zacc1-c1)./abs(zacc1-c1)+tclust(1)]; %the clustering should end at the accumulation point! 
cent_clust2 = [zacc2-tclust'*(zacc2-c2)./abs(zacc2-c2)-tclust(1)];

% cent_clust1 = [zacc1; zacc1-tclust'];
% cent_clust2 = [zacc2; zacc2+tclust'];

cent_clust1 = [zacc1-tclust'*(zacc1-c1)./abs(zacc1-c1)+tclust(1)*(zacc1-c1)./abs(zacc1-c1)]; %the clustering should end at the accumulation point! 
cent_clust2 = [zacc2-tclust'*(zacc2-c2)./abs(zacc2-c2)+tclust(1)*(zacc2-c2)./abs(zacc2-c2)];


%cent_clust1 = [zacc1]; %the clustering should end at the accumulation point! 
%cent_clust2 = [zacc2];

Mclust = Nclust*5; % Possible to pick more clustered collocation points
%Mclust = Nclust*10; %testing
tclust2 = logspace(log10(threshold),0,Mclust)';

%% To assign clustered nodes at the boundary, first find closest parameter value on the two curves
% Objective function: squared distance to zacc1
dist1 = @(t) abs(f1(t)+c1 - zacc1).^2;

% Minimize over t in [0, L]
t_opt1 = fminbnd(dist1, 0, L);

% Closest point on the curve
z_closest1 = f1(t_opt1)+c1;
tclust1 = mod([tclust2; -tclust2]+t_opt1,L);

%repeat for the other obstacle
dist2 = @(t) abs(f2(t)+c2 - zacc2).^2;
t_opt2 = fminbnd(dist2, 0, L);
z_closest2 = f2(t_opt2)+c2;
tclust2 = mod([tclust2; -tclust2]+t_opt2,L);

coll_clust1 = c1 + f1(tclust1);
coll_clust2 = c2 + f2(tclust2);

coll_clust = [coll(1:end/2); coll_clust1; coll(end/2+1:end); coll_clust2];

cent_clust = [centers; cent_clust1; cent_clust2];

%%
figure
plot(cent_clust, 'rx')
hold on
%plot(real(zacc1),imag(zacc1),'b*')
%plot(real(zacc2),imag(zacc2),'b*')
axis equal
plot(coll_clust, 'c+')
title('Enhanced setup')
plot(z_closest1,'k*')
plot(z_closest2,'k*')

%% NOW, solve with the enhanced grid

M = length(coll_clust);
% Get rhs

bvec = [];
Mk = M/P;
start = 1; 
for k = 1:P
   stop = Mk*k;
   rhs_f = @(x) [Up(k,1)-Wp(k)*(imag(x)-imag(q(k))); Up(k,2)+Wp(k)*(real(x-q(k)))];
   rhs_k = rhs_f(coll_clust(start:stop));
   bvec(start:stop) = rhs_k(1:end/2);
   bvec(M+start:M+stop) = rhs_k(end/2+1:end);  
   start = stop+1; 
end
bvec = bvec.';
dir =  rand(2*length(cent_clust1),1);
clust_cent = [cent_clust1; cent_clust2];
[x_c,reserr_c,coefnorm_c] = solve_dense_stokes(centers,coll_clust,clust_cent,dir,svec,bvec,trunc,scalenorms);



%% Compare accuracy

% evaluate accuracy with and without the clustered points. Look at
% different target points than the collocation points

Neval = 1e4;
t_check = linspace(0,L,Neval)'; 

tclust = logspace(-8,0,200)';
tcheck1 = [t_check; mod([tclust; -tclust]+t_opt1,L)];
tcheck2 = [t_check; mod([tclust; -tclust]+t_opt2,L);];

y1 = c1+f1(tcheck1);
y2 = c2+f2(tcheck2);
bndpts_b = [y1; y2]; 

u_p = eval_Stokes_solution(x_b,centers,bndpts_b,[],[],svec,N1,2);


if svec(5)
    dir = [ones(length(clust_cent),1);1i*ones(length(clust_cent),1)];
end
u_c = eval_Stokes_solution(x_c,centers,bndpts_b,clust_cent,dir,svec,N1,2);

M = length(bndpts_b);
% Get rhs

u_ref = [];
Mk = M/P;
start = 1; 
for k = 1:P
   stop = Mk*k;
   rhs_f = @(x) [Up(k,1)-Wp(k)*(imag(x)-imag(q(k))); Up(k,2)+Wp(k)*(real(x-q(k)))];
   rhs_k = rhs_f(bndpts_b(start:stop));
   u_ref(start:stop) = rhs_k(1:end/2);
   u_ref(M+start:M+stop) = rhs_k(end/2+1:end);  
   start = stop+1; 
end
u_ref = u_ref.';


diff_vec_p = u_p-u_ref;
diff_vec_p = diff_vec_p(1:end/2)+1i*diff_vec_p(end/2+1:end);
err_ext_p = max(abs(diff_vec_p));

diff_vec_c = u_c-u_ref;
diff_vec_c = diff_vec_c(1:end/2)+1i*diff_vec_c(end/2+1:end);
err_ext_c = max(abs(diff_vec_c));

fprintf('Without clustering: Residual (collocation points) %1.3e, (new points) %1.3e and coefficient norm %1.3e\n', reserr_p,err_ext_p,coefnorm_p)
fprintf('With clustering: Residual (collocation points) %1.3e, (new points) %1.3e and coefficient norm %1.3e\n', reserr_c,err_ext_c,coefnorm_c)

%% Visualise source strengths -- works only without dipoles


alignfigs;




