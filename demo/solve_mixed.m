function [U,Fres,lambda_xtot,lambda_ytot, lambda_image, rin,rimage,nimage,relres] = solve_mixed(Pw,Pp,q,F,T,image,getres,visualise)
%%solve_mixed solves a mixed resistance and mobility problem for circular
%%objects in 2D Stokes flow (walls and particles).
%
% Syntax: [U,Fres,lambda_tot,rin,rimage,nimage,relres] = solve_mixed(Pw,Pp,q,F,T,getres,image,visualise)
%
% Inputs:
%   Pw    - Number of walls, assumed to be the first Pw objects
%   Pp    - Number of particles, assumeed to be the last Pp objects
%   q     - column vector of complex values coordinates for the objects of size P = Pw+Pp
%   F     - matrix of size Pp x 2 representing the x y coordinates of the of
%           the net forces acting on the particles
%   T     - column vector of size Pp representing the torque acting on the Pp
%           particles
%   image - logical flag to determine if images are in use for close to touching
%         objects and thereby providing an accurate solution. A less accurate
%         solution is determined if image = 0.
%   getres  - logical flag to return relative residual in surface targets
%           different from the collocation points
%   visualise - flag to plot the particle geometry
%
%Outputs:
%   U      - Column vector of computed rigid body motion of size 3Pp x 1
%   Fres   - Column vector of computed net forces and torques of size 3Pw x 1
%   lambda_tot - Column vector of the total source strenghts in all sources,
%         ordered as [lambda_Stokeslets_x, lambda_Stokeslets_y, lambda_image_x,
%         lambda_image_y]
%   rin    - vector of complex valued inner proxy sources on all objects,
%         ordered object by object
%   rimage - struct containing locations for images sources for pair {i,j}
%   nimage - struct containing "directions" (generated at random) for image stresslets for 
%         pair {i,j}
%   relres  - relative residual in points different from the collocation
%           points
%
%Example:
%    See "test_mixed_problem.m"
%
%Notes:
%    Only suitable for small systems, as a dense SVD is used for all
%    objects. Iterative solution would require 2-body preconditioning with
%    different types of 1-body basis functions for resistance and mobility.
%
%Anna Broms April 14, 2025


if nargin==0, test_mixed; return; end

%% Checks

assert(size(T,1)==Pp,'Wrong size of torque vector')
assert(size(F,1)==Pp,'Wrong size of force vector')
assert(size(F,2)==2,'Wrong size of force vector, should contain x y coordinates')
assert(length(q)==Pp+Pw,'Wrong number of center coordinates')

%% PARAMS
P = length(q);

opt = get2Dparams();
opt.image = image;

%Params to set Rp and grid resolution
if image
    Nc = 60; %has been 60 all the time before
   % Nc = 100;
   % Nc = 200;

    tol = 1e-12;
else
    Nc = 200;
    Nc = 1000;
    Nc = 100;

    % more crude
    Nc = 20;
    tol = 1e-6;

    %crude
%     Nc = 50;
%     tol = 1e-12;
    % tol here should probably be larger to allow for smaller coefficients, which in turn will improve on the GMRES convergence
end

sep = (1/Nc)*log(1/tol);
Rp = max([1-sep,0.01]); %radius of proxy surface

%solve with Stresslets +Doublets at image points
s = [0 0 1 1]; % Source types at image points hardcoded
opt.N_c = Nc; 
opt.Rp_c = Rp;
opt.s = s;
opt.pc = 0; 
%% Get and visualise grid of sources and collocation points

%create grids
[rout,rin,rimage,nimage,pair_points] = get2DImageGrid(q,ones(size(q)),opt);

if visualise
    figure(4)
    clf
    %Visualise particles
    plot(real(rin),imag(rin),'b.');
    hold on
    plot(real(rout),imag(rout),'k.');
    plot(real(rimage),imag(rimage),'b.');
    quiver(real(rimage),imag(rimage),real(nimage),imag(nimage))
    axis equal
end

%% Prepare RHS
%get evaluation of lambda0
F = [zeros(Pw,2); F];
T = [zeros(Pw,1); T];
[u,lambda_x,lambda_y] = getRecompletionFlow(rin,rout,q,F,T);

%% Build system matrix for the mixed mobility / resistance problem
solve = 1;
%Determine L first.
[Ntot,Kin,L] = buildMixedSystem(q,rin,rout,rimage,nimage, Nc,Pw,Pp,pair_points,s,solve);

%% Solve densely via SVD
tol = 1e-14;
visualise = 0;
[Y,B]  = getPseudoFactors(Ntot,tol,visualise);

temp = B'*u;
lambda = Y*temp;

%% Post-process to determine velocities on the particles and forces/torques on the walls
x_all = P*Nc;
for k = 1:Pp
    U((k-1)*3+1:3*k) = -Kin'*[lambda(Nc*Pw+(k-1)*Nc+1:Nc*Pw+k*Nc); lambda(x_all+Nc*Pw+1+(k-1)*Nc:x_all+Nc*Pw+k*Nc)];
end
for k = 1:Pw
    Fres((k-1)*3+1:3*k) = Kin'*[lambda((k-1)*Nc+1:k*Nc); lambda(x_all+1+(k-1)*Nc:x_all+k*Nc)];
end

%F1 = Kin'*[lambda(1:Nc); lambda(x_all+1:x_all+Nc)];
%F2 = Kin'*[lambda(Nc+1:2*Nc); lambda(x_all+1+Nc:x_all+2*Nc)];
%Should be smaller for smaller distance? But, Stokes paradox?

%check residual in SAME set of points as the collocation points.
str = sprintf('Rel residual in sol, at colloc points %1.2e, Accurate scheme? %u',norm(Ntot*lambda-u),image);
disp(str)

%Determine total source strengths (undo the splits necessary for mobility).
%
% source strengths on walls
lambda_xtot = lambda(1:Nc*Pw);
lambda_ytot = lambda(P*Nc+1:P*Nc+Nc*Pw);

% get image sources
lambda_image = lambda(P*Nc*2+1:end); 

% project source strengths on particles
for k = 1:Pp
    lambda_k = (eye(size(L))-L)*[lambda(Nc*Pw+1+(k-1)*Nc:Nc*Pw+k*Nc); lambda(x_all+Pw*Nc+1+(k-1)*Nc:x_all+Pw*Nc+k*Nc)];
    lambda_xtot = [lambda_xtot; lambda_k(1:end/2)+lambda_x(Nc*Pw+(k-1)*Nc+1:Nc*Pw+k*Nc)];
    lambda_ytot = [lambda_ytot; lambda_k(end/2+1:end)+lambda_y(Nc*Pw+(k-1)*Nc+1:Nc*Pw+k*Nc)];
end

lambda_tot = [lambda_xtot; lambda_ytot; lambda_image];


%% Also, check in new points, using direct summation
if getres
    np = 1700; %evaluate residual in a large number of points per body
    rcheck = [];
    tt = linspace(0,2*pi,np+1);
    tt = tt(1:end-1)';
    Klarge = getKmat2D(cos(tt)+1i*sin(tt),0);
    for k = 1:P
        rcheck = [rcheck; q(k)+cos(tt)+1i*sin(tt)];
        if k > Pw
            uk = Klarge*U((k-Pw-1)*3+1:3*(k-Pw))';
            u_RBM((k-1)*np+1:k*np) = uk(1:end/2); 
            u_RBM(P*np+(k-1)*np+1:P*np+k*np) = uk(end/2+1:end);
        end
    end
    u_RBM = u_RBM';
    stress_x = lambda_image(1:end/4);
    stress_y = lambda_image(end/4+1:end/2);
    pot_x = lambda_image(end/2+1:3*end/4);
    pot_y = lambda_image(3*end/4+1:end); 
    u = getVelocityField(rin,rcheck,lambda_xtot,lambda_ytot,rimage,nimage,[],stress_x,stress_y,pot_x,pot_y);
    diff_vec = u-u_RBM;
    max_val = max(abs(u_RBM(1:end/2)+1i*u_RBM(end/2+1:end)));
    relres = max(abs(diff_vec(1:end/2)+1i*diff_vec(end/2+1:end)))/max_val;

else
    relres = [];
end




end

function [Ntot,Kin,L] = buildMixedSystem(q,rin,rout,rim,nim,Nc,Pw,Pp,pair_points,s,solve)
%Builds system matrix for a mixed resistance and mobility problem. Nc is
%the number of coarse source points on the proxy boundary of each particle.
% pair_points contains in its second column and row j information on the number of
% collocation points for particle j. s is a vector that sets the
% singularity types along image lines.

%Image enhancement?
if size(rim)
    image = 1;
else
    image = 0;
end

%In the system matrix, get Stokeslets without projection for every particle
%where a velocity is prescribed (all walls)
Nio = singleLayer(rin(1:Pw*Nc),rout,1);
%The first objects are fixed walls and the last objects are assumed to be moving particles

%Create projection for the particle(s) for which a force is provided (the
%moving bodys)
rin_k = rin(1:Nc);
Kin = getKmat2D(rin_k,q(1));
L = Kin*((Kin'*Kin)\Kin'); %Projection matrix, same for everyone as it depends only on the proxy sources

M = sum(pair_points(1:Pw,2));
Mtot = sum(pair_points(:,2));

%Loop over all particles to close system (express unknown RBM velocities in
%terms of the unknown sources). 
Lr_all = zeros(2*sum(pair_points(:,2)),2*Nc*(Pp));
for k = 1:Pp
    rout_k = rout(M+1:M+pair_points(k+Pw,2)); %get collocation points on particles!
    Kout = getKmat2D(rout_k,q(k+Pw)); 
    Lr = Kout*Kin'; %particle dependent!

    Lr_all(M+1:M+pair_points(k+Pw,2),Nc*(k-1)+1:Nc*k) = Lr(1:pair_points(k+Pw,2),1:end/2);
    Lr_all(M+1:M+pair_points(k+Pw,2),Nc*Pp+Nc*(k-1)+1:Nc*Pp+Nc*k) = Lr(1:pair_points(k+Pw,2),end/2+1:end);

    Lr_all(Mtot+M+1:Mtot+M+pair_points(k+Pw,2),Nc*(k-1)+1:Nc*k) = Lr(pair_points(k+Pw,2)+1:end,1:end/2);
    Lr_all(Mtot+M+1:Mtot+M+pair_points(k+Pw,2),Nc*Pp+Nc*(k-1)+1:Nc*Pp+Nc*k) = Lr(pair_points(k+Pw,2)+1:end,end/2+1:end);
    M = sum(pair_points(1:Pw+k,2));
end

start_colloc = 0;

P = Pp+Pw;

% Project so that Stokeslet sources on the particle proxy surfaces don't
% contribute to particle net force and torque
for k = 1:P % Loop over targets on everyone
    M_k = pair_points(k,2);
    rout_k = rout(start_colloc+1:start_colloc+M_k);
    %Loop over particles
    for i = 1:Pp
        Npk = singleLayer(rin((Pw+(i-1))*Nc+1:(Pw+i)*Nc),rout_k,1);
        A = Npk*(eye(size(L))-L);

        last_col_blocks(start_colloc+1:start_colloc+M_k,(i-1)*Nc+1:i*Nc) = A(1:end/2,1:end/2);
        last_col_blocks(start_colloc+1:start_colloc+M_k,Pp*Nc+(i-1)*Nc+1:Pp*Nc+i*Nc) = A(1:end/2,end/2+1:end);
        last_col_blocks(start_colloc+1+sum(pair_points(:,2)):start_colloc+sum(pair_points(:,2))+M_k,(i-1)*Nc+1:i*Nc) = A(end/2+1:end,1:end/2);
        last_col_blocks(start_colloc+1+sum(pair_points(:,2)):start_colloc+sum(pair_points(:,2))+M_k,Pp*Nc+(i-1)*Nc+1:Pp*Nc+i*Nc) = A(end/2+1:end,end/2+1:end);
    end
    start_colloc = start_colloc+M_k;
end

if solve
    last_col_blocks = last_col_blocks+Lr_all;
end
Ntot = [Nio(:,1:end/2) last_col_blocks(:,1:end/2) Nio(:,end/2+1:end) last_col_blocks(:,end/2+1:end)];

%%


if image
    Nimage = getImageKernels2D(rim,nim,rout,1,s);
else
    Nimage = [];
end
Ntot = [Ntot Nimage];

end

function test_mixed
    close all; 
    Pw = 2;
    Pp = 1;

    %Set force and torque on paritlce
    F = [1 0];
    T = 0;

    %positions for the objects (not close)
    a = 4;
    delta = 0.001;
    y = [0 0 -a; 2+delta -2-delta 0];
    q = y(1,:)+1i*y(2,:);
    image = 1; %will use a finer grid
    [Uacc,Facc,lambda_acc,rin_acc,rimage,nimage,s] = solve_mixed(Pw,Pp,q,F,T,image,1);
    image = 0; %coarser grid
    [Ucoa,Fcoa,lambda_coa,rin_coa,rimage,nimage,s] = solve_mixed(Pw,Pp,q,F,T,image,1);

    disp('Relative difference in computed velocities, dilute setting')
    norm(Ucoa-Uacc)/norm(Uacc)
    %The objects are well-separated and the accuracy is high also with a
    %coarse grid!

    %Try to compute a solution also with the three objects in a line (a
    %much harder problem!)
    a = 0;
    delta = 0.001;
    y = [0 0 -a; 2+delta -2-delta 0];
    q = y(1,:)+1i*y(2,:);
    image = 1; 
    [Uacc,Facc,lambda_acc,rin_acc,rimage,nimage,s] = solve_mixed(Pw,Pp,q,F,T,image,1);
    image = 0; 
    [Ucoa,Fcoa,lambda_coa,rin_coa,rimage,nimage,s] = solve_mixed(Pw,Pp,q,F,T,image,1);
    disp('Relative difference in computed velocities, particles in a row')
    norm(Ucoa-Uacc)/norm(Uacc)



end
