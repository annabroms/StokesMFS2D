function [err_extp,err_extl,reserr_p,reserr_l,coefnorm_p,coefnorm_l] = june21_smooth_stokes(svec,N,plotdomain)

if nargin<1
    close all
    N = 150; 
    svec = [0 0 0 1 1];
    svec = [1 1 0 0 0];
    %svec = [1 0 0 0 0];
    plotdomain = 0;
    
end


%% Setup a problem with an ellipse

% This is an experiment to test a Stokes problem, similar to the june14
% script for Laplace and the Helmholtz 20 June script 

fprintf('------Experiment with 2D Stokes ------\n\n')

a = 2.0;
b = 0.2;
%b = 1; 
c = sqrt(a^2-b^2);   % foci of ellipse are at (+/-c, 0)

osf = 10;
%osf = 30; 
regtol = 1e-13;
scalecols = 1;


%N = 400; 
M = osf*N;
r_proxy = c/a;
tau = 0.4*b;


plotfield = 0;



t = linspace(0, 2*pi, 401);
z1 = a*cos(t)+b*1i*sin(t);

t_proxy = linspace(0, 2*pi, N+1);
t_proxy(end) = [];
%centers_proxy = (a*cos(t_proxy)+b*1i*sin(t_proxy))*r_proxy;
centers_proxy = (a*cos(t_proxy+1i*tau)+b*1i*sin(t_proxy+1i*tau)); %better choice of disc
centers_proxy = centers_proxy.';

centers_line = chebpts(round(N/2))' * r_proxy * a;
% centers_line = linspace(-r_proxy*a, r_proxy*a, round(N/2)); worse choice
centers_line = centers_line.';

tau = linspace(0, 2*pi, M+1);
tau(end) = [];
coll = a*cos(tau)+b*1i*sin(tau);
coll = coll.';

%% build rhs 
U = [1 2]; %translational velocity on the entire square
%perhaps do rotation and translation?
c = 0;  
W = 2; 
rhs_f = @(x) [U(1,1)-W*(imag(x)-imag(c)); U(1,2)+W*(real(x-c))]; %quickly varying data for a slender ellipse
%rhs_f = @(z) [real(z).^2; -real(z).^2]; %not diverengce free! 
%rhs_f = @(z) [real(z).^2; -2*real(z).*imag(z)]; %divergence free, needed for the interior problem.
bvec = rhs_f(coll); 

%Get check points
tau = linspace(0, 2*pi, 1e4+1);
tau(end) = [];
bndpts_b = a*cos(tau)+b*1i*sin(tau);
bndpts_b = bndpts_b.';
u_ref = rhs_f(bndpts_b); 

%%

if plotdomain
    [r, pol, res, zer, zj, fj, wj, errvec] = aaa(conj(z1), z1,'tol',1e-6);
    plot(real(z1), imag(z1), 'k-')
    axis equal
    hold on
    plot(real(centers_proxy), imag(centers_proxy), 'rx')
    plot(real(centers_line), imag(centers_line), 'bx')
    plot(real(pol),imag(pol),'ko') %look at singularities of Schwarz function
    xlim([-3,3])
    ylim([-3,3])

    figure()
    scatter(real(coll),imag(coll),5,bvec(1:end/2),'filled')
    axis equal
    colorbar
    title('Boundary data')
end

%% Solve using the proxy surface
%dir = 1i*ones(length(centers_proxy),1);
dir = rand(length(centers_proxy),1)+1i*rand(length(centers_proxy),1);
[x_p,reserr_p,coefnorm_p] = solve_dense_stokes([],coll,centers_proxy,dir,svec,bvec,regtol,scalecols);

if svec(5)
    dir = [ones(length(centers_proxy),1);1i*ones(length(centers_proxy),1)];
end
u = eval_Stokes_solution(x_p,[],bndpts_b,centers_proxy,dir,svec,0);
diff_vec = u-u_ref;
diff_vec = diff_vec(1:end/2)+1i*diff_vec(end/2+1:end);
err_extp = max(abs(diff_vec));

fprintf('Proxy surface: residual %1.3e and coefficient norm %1.3e\n\n', reserr_p, coefnorm_p)


%% Solve using an interval of centers

dir = 1i*ones(length(centers_line),1);
dir = rand(length(centers_line),1)+1i*rand(length(centers_line),1);
[x_l,reserr_l,coefnorm_l] = solve_dense_stokes([],coll,centers_line,dir,svec,bvec,regtol,scalecols);

if svec(5)
    dir = [ones(length(centers_line),1);1i*ones(length(centers_line),1)];
end
u = eval_Stokes_solution(x_l,[],bndpts_b,centers_line,dir,svec,0);
diff_vec = u-u_ref;
diff_vec = diff_vec(1:end/2)+1i*diff_vec(end/2+1:end);
err_extl = max(abs(diff_vec));

fprintf('Line sources: residual %1.3e and coefficient norm %1.3e\n\n', reserr_l, coefnorm_l)


%% Plot the field

if plotfield
    Xrange = -4:0.05:4;
    Yrange = -4:0.05:4;
    [X,Y] = meshgrid(Xrange, Yrange);
    
end


end