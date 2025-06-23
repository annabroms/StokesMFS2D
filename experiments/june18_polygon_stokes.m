function [err_ext,reserr_ext,coefnorm_ext] = june18_polygon_stokes(svec,plotdomain)
% Investigates Stokes MFS for a domain with corners. The input svec
% determines the set of singularities clustered towards the corners. The
% interior problem is much easier to resolve than the exterior one, and no
% clustered nodes are needed in the interior setting. The non-convex
% exterior problem is particularly challenging. In the script
% "june19_polygon_stokes_singularities.m", we explore suitable
% singularities at the clustered nodes. To visualise example, run this
% function without arguments. 
%
% To do: play with row weighting at the solve stage 

if nargin<1
    close all
    svec = [1 1 0 1 1]; %S R Tr, D, T. [1 1 0 1 1] works for the non-convex polygon exterior/interior problem
                        %Tr refers to Stresslet with random "normal"
                        %direction 
   % svec = [0 0 0 1 1];                 
    plotdomain = 1;
elseif nargin<2
    plotdomain = 0; 
end

regtol = 1e-13;
%regtol = 1e-10; %truncation level of SVD

solve_int = 0; %solve interior problem? 

%% Definition of the polygon

P1 = [1+2i 1/2+1i/2 1+2/3+2i/3 2+2i 3/2+3i/2];
P2 = [1+2i 1/2+1i/2 1+2/3+2i/3 2+2i];

%P2 = [2i 0 2 2+2i]; %orientation is important here for orientation of normals 
square = 0; 

P = P1;  % for a non-convex polygon
%P = P2;  % for a convex polygon

n_sides = length(P);

% compute tangential and normal directions
P_per = [P P(1)];
diffP = diff(P_per);
tangents = diffP ./ abs(diffP);
normals = exp(1i*(angle(tangents)-pi/2));
normals_per = [normals(end) normals];


% the outward normals at the vertices
vert_normals = 1/2 * (normals_per(1:n_sides) + normals_per(2:n_sides+1));
vert_normals = vert_normals ./ abs(vert_normals);
vert_tangents = exp(1i*(angle(vert_normals)-pi/2));

%% Parameters of the problem
%n_smooth = 200;
n_smooth = 150; 
n_smooth = 100; 
n_clust = 150; %150 takes a long time... 
%n_clust = 0; %the interior problem is very easy to solve. How come?
%n_clust = 10; 


%n_clust = 100; 
A_clust = 0.6; %original setting
sigma = 4; %Will impact accuracy!
%sigma = 2; 

osf = 3;    % oversampling factor

if square
    r_proxy = 0.6;
    r_proxy = 0.7;
else
    r_proxy = 0.25;
    r_proxy = 0.35;
    r_proxy = 0.15;
end

%% Computation of the centers

% the clustering centers
s = sqrt(1:n_clust) - sqrt(n_clust);
d = A_clust * exp(sigma*s);
%d = d(d>1e-15);
d = d(d>1e-14);
%d = d(d>1e-12);

%For n_clust = 150, d = d(d>1e-15) gives 
%large residual in new points for the interior
%problem! The residual at the collocation points is
%still small. Everything works well for smaller n_clust. If I change d =
% d(d>1e-15); to d = d(d>1e-14); all seems fine.

centers_clust = [];
dirs = [];

for k = 1:n_sides
    c = P(k) - vert_normals(k)*d;
    centers_clust = [centers_clust c];
    dirs = [dirs vert_tangents(k)*ones(1,length(c))];       
end

centers_clust_int = [];
for k = 1:n_sides
    c = P(k) + vert_normals(k)*d;
    centers_clust_int = [centers_clust_int c];
end

% the smooth centers -- will impact acuracy!
theta = linspace(0, 2*pi, n_smooth+1)';
theta(end) = [];

c = 1.2+1.2i; % proxy center
c = 1.25+1.1i;
centers_smooth = c + r_proxy * exp(1i*theta);
centers_smooth_int = c + 2 * exp(1i*theta);

% works well for the non-convex polygon
centers_smooth = c + r_proxy * exp(1i*theta);
centers_smooth_int = c + 2 * exp(1i*theta);

% centers_smooth = mean(P) + r_proxy * exp(1i*theta);
% centers_smooth_int = mean(P) + 2 * exp(1i*theta);


if square
    centers_smooth = 1+1i+ r_proxy * exp(1i*theta);
    centers_smooth_int = 1+1i+ 3* exp(1i*theta);
end



centers = [centers_clust centers_smooth.'].';
centers_int = [centers_clust_int centers_smooth_int.'].';

%% Computation of boundary points
pts1 = clustered_points(round(n_smooth/n_sides), n_clust, sigma);
pts = oversample_points(pts1, osf);
pts = [-1 pts 1]; %add corner points 

bndpts = [];
weights = []; 
for k = 1:n_sides
    pk = (P_per(k)+P_per(k+1))/2 + pts/2 * abs(P_per(k+1)-P(k)) * tangents(k);
    w = min(abs(pk - P_per(k)), abs(pk - P_per(k+1)));
    bndpts = [bndpts pk];
    w(1) = 1; %w(end) = 1; % modify weight at corners
    weights = [weights w]; 
end
bndpts = bndpts.';
weights = weights';
%weights = weights.^2; % a la Trefehten's lightning solver
weights = ones(size(weights)); %turn off row weighting. How to properly deal with the corners?

%% Plot the domain

if plotdomain
    figure
    plot(bndpts, 'k-')
    hold on
    plot(centers, 'ro')
    plot(centers_int, 'gx')
    hold off
    axis equal
end

%%

scalecols = 1;
dir = (1+1i)*ones(length(centers_clust),1);

%build rhs 
U = [1 2]; %translational velocity on the entire square
%perhaps do rotation and translation?
c = mean(P); %rotation center -- not sure this is sensible...  
W = 2; 
rhs_f = @(x) [U(1,1)-W*(imag(x)-imag(c)); U(1,2)+W*(real(x-c))];
%rhs_f = @(z) [real(z).^2; -real(z).^2]; %not diverengce free! 
%rhs_f = @(z) [real(z).^2; -2*real(z).*imag(z)]; %difvergence free, needed for the interior problem.
b = rhs_f(bndpts); 

% Two alternatives here are:
% - at the clustered centers, we use both the SLP kernel and the DLP kernel
% in the normal direction

centers_clust = centers_clust.'; 
centers_clust_int = centers_clust_int.';

%% A different set of boundary points

pts_b = oversample_points(pts1, osf+1);
pts_b = [-1 pts_b 1]; %add corner points 

bndpts_b = [];
for k = 1:n_sides
    pk = (P_per(k)+P_per(k+1))/2 + pts_b/2 * abs(P_per(k+1)-P(k)) * tangents(k);
    bndpts_b = [bndpts_b pk];
end

bndpts_b = bndpts_b.';

u_ref = rhs_f(bndpts_b);

%% Solve exterior problem

[x_out,reserr_ext,coefnorm_ext,solres] = solve_dense_stokes(centers_smooth,bndpts,centers_clust,dir,svec,b,regtol,scalecols,weights);

if svec(5)
    dir = [ones(length(centers_clust),1);1i*ones(length(centers_clust),1)];
end
u = eval_Stokes_solution(x_out,centers_smooth,bndpts_b,centers_clust,dir,svec,n_smooth);
diff_vec = u-u_ref;
diff_vec_c = diff_vec(1:end/2)+1i*diff_vec(end/2+1:end);
err_ext = max(abs(diff_vec_c));


if plotdomain
    %Look at solution
    % figure()
    % semilogy(real(bndpts_b),abs(u(1:end/2)))
    % hold on
    % semilogy(real(bndpts_b),abs(u_ref(1:end/2)))
    % semilogy(real(bndpts_b),abs(u(end/2+1:end)))
    % semilogy(real(bndpts_b),abs(u_ref(end/2+1:end)))

    figure()
    subplot(1,2,1)
    scatter(real(bndpts),imag(bndpts),5,log10(abs(solres(1:end/2))),'filled')
    colorbar
    subplot(1,2,2)
    scatter(real(bndpts),imag(bndpts),5,log10(abs(solres(end/2+1:end))),'filled')
    colorbar
    sgtitle('Residual at collocation nodes -> zoom in on corners!')

    figure()
    subplot(1,2,1)
    scatter(real(bndpts_b),imag(bndpts_b),5,log10(abs(diff_vec(1:end/2))),'filled')
    title('ux')
    colorbar
    subplot(1,2,2)
    scatter(real(bndpts_b),imag(bndpts_b),5,log10(abs(diff_vec(end/2+1:end))),'filled')
    colorbar
    title('uy')
    sgtitle('Residual at new nodes')

    figure()
    subplot(2,2,1)
    semilogy(real(bndpts_b),abs(diff_vec(1:end/2)),'.')
    xlabel('Re(z)')
    title('ux')
    subplot(2,2,2)
    semilogy(real(bndpts_b),abs(diff_vec(end/2+1:end)),'.')
    xlabel('Re(z)')
    title('uy')
    sgtitle('Residual at new nodes')
    subplot(2,2,3)
    semilogy(imag(bndpts_b),abs(diff_vec(1:end/2)),'.')
    xlabel('Im(z)')
    title('ux')
    subplot(2,2,4)
    semilogy(imag(bndpts_b),abs(diff_vec(end/2+1:end)),'.')
    xlabel('Im(z)')
    title('uy')
    sgtitle('Residual at new nodes')

end


%% Solve interior problem

if solve_int
    [x,reserr_int,coefnorm_int] = solve_dense_stokes(centers_smooth_int,bndpts,centers_clust_int,dir,svec,b,regtol,scalecols);
    
    u = eval_Stokes_solution(x,centers_smooth_int,bndpts_b,centers_clust_int,dir,svec,n_smooth);
    
    diff_vec = u-u_ref;
    diff_vec_c = diff_vec(1:end/2)+1i*diff_vec(end/2+1:end);
    err_int = max(abs(diff_vec));
    
    % figure()
    % semilogy(abs(u))
    % hold on
    % semilogy(abs(u_ref))

    figure()
    subplot(1,2,1)
    scatter(real(bndpts_b),imag(bndpts_b),5,log10(abs(diff_vec(1:end/2))),'filled')
    colorbar
    subplot(1,2,2)
    scatter(real(bndpts_b),imag(bndpts_b),5,log10(abs(diff_vec(end/2+1:end))),'filled')
    colorbar
end


%%
fprintf('Exterior problem: residual %1.3e and coefficient norm %1.3e, residual in new points %1.3e\n', reserr_ext, coefnorm_ext,err_ext)

if solve_int
    fprintf('Interior problem: residual %1.3e and coefficient norm %1.3e, residual in new points %1.3e\n', reserr_int, coefnorm_int,err_int)
end


if plotdomain
    figure
    plot(bndpts, 'k-')
    hold on
    scatter(real(centers_smooth),imag(centers_smooth),20,log10(abs(x_out(1:length(centers_smooth)))),'filled')
    scatter(real(centers_clust),imag(centers_clust),20,log10(abs(x_out(end-length(centers_clust):end-1))),"filled")
    axis equal
    colorbar
    title('Source strength magnitude log10, basic')
    hold off

    figure()
    Npts = 5000; %number of evaluation points
    pts = -1-1i+5*(rand(Npts,1)+1i*rand(Npts,1));
    [in,on] = inpolygon(real(pts),imag(pts),real(P),imag(P));
    
    lim = 4; %set clim
  
    ext_points = [pts(~in); bndpts_b];
    u = eval_Stokes_solution(x_out,centers_smooth,ext_points,centers_clust,dir,svec,n_smooth);
    subplot(1,2,1)
    surfir(real(ext_points),imag(ext_points),u(1:end/2));
    hold on
    colorbar
    clim([-lim,lim])
    fill3(real(P),imag(P), 20*ones(size(P)), [1 1 1]); %'w', 'FaceAlpha', 100, 'EdgeColor', 'w');
    title('x-velocity')
    view(0,90)
    xlim([-0.5,3])
    ylim([-0.5 3.5])
    axis square

    subplot(1,2,2)
    surfir(real(ext_points),imag(ext_points),u(end/2+1:end));
    hold on
    colorbar
    clim([-lim,lim])
    fill3(real(P),imag(P), 20*ones(size(P)), [1 1 1]);
    title('y-velocity')
    view(0,90)
    xlim([-0.5,3])
    ylim([-0.5 3.5])
    axis square 
    sgtitle('Exterior problem')

    %visualise evaluation points
    figure()
    plot(real(ext_points),imag(ext_points),'*')
    title('Evaluation points')

    %% Also visualise the interior solution
    if solve_int
        int_points = [pts(in); bndpts_b];
        u = eval_Stokes_solution(x,centers_smooth_int,int_points,centers_clust_int,dir,svec,n_smooth);  
        subplot(1,2,1)
        surfir(real(int_points),imag(int_points),u(1:end/2));
        hold on
        colorbar
        clim([-lim,lim])
        title('x-velocity')
        view(0,90)
        xlim([0,2.5])
        ylim([0 2.5])
        axis square
    
        subplot(1,2,2)
        surfir(real(int_points),imag(int_points),u(end/2+1:end));
        hold on
        colorbar
        clim([-lim,lim])
        title('y-velocity')
        view(0,90)
        xlim([0,2.5])
        ylim([0 2.5])
        axis square
        sgtitle('Interior problem')
    
        
    end
    alignfigs();

end

end


