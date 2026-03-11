function [UW,lambda_c,it,gmres_tol, rel_res,abs_res] = solve_mob_peanut_enhanced(q,F,T,delta_pair,N_peanut,visualise,gmres_tol,debug,surface_error_mode,gmres_verbose)
%SOLVE_MOB_PEANUT_ENHANCED Solve a 2D Stokes mobility problem with
%peanut-compressed pair corrections. The enhanced discretisation for each
%pair consists of Stokeslets only. 
%
% Syntax:
%   [UW,lambda,it,gmres_tol,rel_res,abs_res] = solve_mob_peanut_enhanced( ...
%       q,F,T,delta_pair,N_peanut,visualise,gmres_tol,debug)
%
% Inputs:
%   q          - Complex particle centers (P x 1).
%   F          - Net forces per particle (P x 2), columns are x/y.
%   T          - Net torques per particle (P x 1).
%   delta_pair - Near-pair threshold; pairs closer than this use local fine correction.
%   N_peanut   - Number of peanut check points used in pair compression.
%   visualise  - Plot diagnostics if true.
%   gmres_tol  - Optional GMRES tolerance (default 1e-10)
%   debug      - Optional logical flag: build/draw dense matrix generated 
%                by matvec and its
%                eigenvalues (default false)
%   surface_error_mode - Optional boundary-error plot mode: 'abs' (default)
%                or 'rel'
%
% Outputs:
%   UW         - Rigid-body motion vector [Ux1; Uy1; W1; ...; UxP; UyP; WP].
%   lambda_c   - Stacked coarse source strengths [lambda_x; lambda_y]
%                based on compressed grid.
%   it         - GMRES iterations.
%   gmres_tol  - GMRES tolerance used.
%   rel_res    - Maximum relative boundary residual on an independent check grid.
%   abs_res    - Maximum absolute boundary residual on an independent check grid.
%
% Notes:
%   - Call with no inputs to run the built-in test.
%   - All Stokeslet sources contribute to the rigid body velocities
%     computed in UW.
%
% Anna Broms, Feb 23, 2026

if nargin==0, test_solve_mob; 
    return; end

if nargin < 6 || isempty(visualise), visualise = 0; end
if nargin < 7 || isempty(gmres_tol), gmres_tol = 1e-10; end
if nargin < 8 || isempty(debug), debug = false; end
if nargin < 9 || isempty(surface_error_mode), surface_error_mode = 'rel'; end
if nargin < 10 || isempty(gmres_verbose), gmres_verbose = 0; end
surface_error_mode = lower(char(surface_error_mode));
if ~any(strcmp(surface_error_mode, {'abs','rel'}))
    error('surface_error_mode must be ''abs'' or ''rel''.')
end

%% SET PARAMS
%GMRES params
maxit = 800; % max GMRES iterations
solver_name = 'mob_peanut_enhanced';

%Set coarse and fine grid. 
%Play with N_c, N_f, a (a_f). 
N_c = 80;  %Number of Stokeslets per body on coarse grid.
% higher resolution better here for dense systems
%N_c = 150; 

N_f = 150; %Number of Stokeslets per body on coarse grid. Only used for the 2-body BVP
%N_f = N_c; 

%N_c = 250; 
a_c = 1.2; % ratio between number of source and collocation points for coarse proxy grid
a_f = 1.2; % ... for the fine proxy grid

tol_c = 1e-12; %I think this works reasonably
tol_c = 1e-10; %Curve moves closer to the surface -> smaller coeff  
%tol_c = 1e-16; %Curve moves further from surface -> larger coeff. 

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rp_c = max([1-sep_c,0.01]); %radius of proxy surface for coarse grid
Rp_f = max([1-sep_f,0.01]);  % and fine grid
%Rp_f = Rp_c; %testing

%radius for which additional, enhancing nodes, should not be used. Based on
%closed formula from fixed point of reflection.
accstop = (1-Rp_c)^2/Rp_c;  

if nargin < 4 || isempty(delta_pair)
    delta_pair = accstop; %We want to use the pair correction for all gaps smaller than delta_pair. (or accstop). 
    % We may choose delta_pair larger than accstop.
end

opt = get2Dparams(); 
opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c; 
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 
opt.N_peanut = N_peanut; %numer of points on peanut boundary
opt.precomp = 1; %faster if evaluation of one body basis on fine grid is compted only once. 
% %Less storage required.
opt.pc = 1; %prepare grid to do pair corrections
opt.delta_pair = delta_pair; 
opt.bndry_vel = 1; 

opt.cmap = 1; % Use coarse to coarse mapping by multiplying together all factors needed for the fine grid mapping 
opt.Nclust = 100; %points on ellipse segments, for now (not the actual number as some will fall inside the proxy curve)

P = length(q);
opt.P = P; 
opt.gmres_verbose = gmres_verbose;
%% Discretize
%Outer basic grid
tout_c_all = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c_all(1:end-1)';
rbase_out_c = cos(tout_c)+1i*sin(tout_c);
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rp_c*cos(tin)+Rp_c*1i*sin(tin); 

tin_f = linspace(0,2*pi,N_f+1);
tin_f = tin_f(1:end-1)';
rbase_in_f =  Rp_f*cos(tin_f)+Rp_f*1i*sin(tin_f);

[~, ~, ~, rimage_vec, refine,pairs] = getEnhancedGrid(q, opt);

rvec_in_c = [];
rout = [];

for k = 1:P
    rvec_in_c = [rvec_in_c; q(k)+rbase_in_c];
    rout = [rout; rbase_out_c+q(k)];
end

%% Prepare completion flow, one- and two-body factorisations
%get completion flow: evaluation of lambda0
[urhs,~] = getRecompletionFlow(rvec_in_c,rout,q,F,T); 

%Get one-body pseduo inverse blocks -- enough to do this for single body.
rimage_in = []; 
[U,Y,Lc] = getSelfPseudoMobilityStokes(1,q,rbase_in_c,rbase_out_c,rimage_in,[0,ceil(a_c*N_c)]);

%Get pair basis
plot_grid = 0; 
opt.project = true;
opt.pair_basis_debug = plot_grid;
opt.show_counter = true;
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU] = getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,Lc{1});
                              

% TODO: update visualisation:
%Visualise 1-body and pair-basis
%warning('Deactivate opt.precomp');
%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,[],[],N_c, N_f,a_c,a_f,rads)

Lc_pair = getILpair(Lc{1}); %I-L projection with coarse grid for a pair. TODO: replace.

%% Construct check boundaries
%rcheck_dom = 100+100i; %some point far away

% Create new grid points, for which the accuracy of the solution is
% to be evaluated. 
rcheck_b = [];
n_bound = 803;
t = linspace(0,2*pi,n_bound)';
for k = 1:P
    rcheck_b = [rcheck_b; q(k)+cos(t)+1i*sin(t)];
end

%% Repackage to prepare for solve
% Shared input structs for matvec and transform calls.
% Use dot-assignment to keep cell fields intact and avoid struct-array expansion.
geom_solve = struct();
geom_solve.rbase_in_c = rbase_in_c;
geom_solve.rbase_in_f = rbase_in_f;
geom_solve.rvec_in = rvec_in_c;
geom_solve.refine = refine;
geom_solve.rimage_vec = rimage_vec;
geom_solve.opt = opt;
geom_solve.rvec_out = rout;
geom_solve.rcheck = rout;
geom_solve.q = q;
geom_solve.pairs = pairs;

basis_mob = struct();
basis_mob.U = U;
basis_mob.Y = Y;
basis_mob.Lc = Lc{1};
basis_mob.Upf = UB_all;
basis_mob.Ypf = YB_all;
basis_mob.DC_all = UC_all;
basis_mob.YC_all = YC_all;
basis_mob.Cmap = Cmap;
basis_mob.Cmap_FU = Cmap_FU; 
basis_mob.Lc_pair = Lc_pair;

geom_check = geom_solve;
geom_check.rcheck = rcheck_b;

%% Solve system

% Debug mode: build the matrix to check it out
if debug
    x = zeros(2*length(rout),1);
    tic
    for k = 1:2*length(rout)
        k
        x(:) = 0; 
        x(k) = 1; 
        uu = matvec_mob_peanut_enhanced(x,geom_solve,basis_mob);
        CC(:,k) = uu;
    end
    toc
    figure(14);
    clf; 
    imagesc(log10(abs(CC)))
    colorbar
    title([solver_name ': log_{10} |CC|'],'interpreter','none')
    skeel(CC)

    [V,D] = eig(CC);
    D = diag(D); 
    figure()
    plot(real(D),imag(D),'+')
    xlabel('Re \lambda')
    ylabel('Im \lambda')
    title([solver_name ': eigenvalues of CC'],'interpreter','latex')
end

% Solve
[tau,it,resvec,real_res] = helsing_gmres(@(x) matvec_mob_peanut_enhanced(x,geom_solve,basis_mob),...
    urhs,2*size(rout,1),maxit,gmres_tol,opt,rout);

figure()
semilogy(resvec); 
title('GMRES convergence with peanut compression, mobility', 'interpreter','latex')


if visualise
    %check residual
    restot = (matvec_mob_peanut_enhanced(tau,geom_solve,basis_mob)-urhs)./urhs;
    figure()
    semilogy(abs(restot))
    title([solver_name ': Rel res at colloc points'],'interpreter','latex')
end
% hold on
% semilogy(resvec2); 

%% COMPUTE Rigid body motion
%And evaluate residual in new points rcheck_b

% Recover coarse and fine sources from data on the boundary
[lam_c_x, lam_c_nonpx,lam_self_x, lam_f_x,lam_c_y, ...
    lam_c_nonpy,lam_self_y,lam_f_y,~,rimage_k] = ...
    transform_mob_peanut_stokes(tau,geom_check,basis_mob);
lambda_c = [lam_c_x; lam_c_y];

%%% Get rigid body motion. 

%First due to all coarse sources
Kc = getKmat2D(rbase_in_c,0);
UW= zeros(3*P,1); 
warning('are the right sources used here?')
for k= 1:P
    UW((k-1)*3+1:3*k) = -Kc'*[lam_c_nonpx((k-1)*N_c+1:k*N_c); lam_c_nonpy((k-1)*N_c+1:k*N_c)];
end

if opt.cmap
    for pair_it = 1:size(pairs,1)
        i = pairs(pair_it,1);
        p2 = pairs(pair_it,2);

        coarse_i = (i-1)*N_c+1:i*N_c;
        coarse_p2 = (p2-1)*N_c+1:p2*N_c;

        warning('check sign in formulae')

        %Here the projected sources are used, as (I-L) is not yet
        %applied for Cmap_FU.
        rhs_pair = [lam_self_x(coarse_i); lam_self_x(coarse_p2); ...
                    lam_self_y(coarse_i); lam_self_y(coarse_p2)];

        % Determine rigid body motion for the pair, using ansatz
        pair_vel = Cmap_FU{i,p2}*rhs_pair;
        UW((i-1)*3+1:3*i) = UW((i-1)*3+1:3*i)+ pair_vel(1:3); 
        UW((p2-1)*3+1:3*p2) = UW((p2-1)*3+1:3*p2)+ pair_vel(4:6);

    end

else
    % Then, due to all pair sources (fine-body + fine-image).
    has_neigh = sort(unique(pairs(:)));
    for i = 1:length(has_neigh)
        k = has_neigh(i);
        rsrc_k = [rbase_in_f+q(k); rimage_k{k}];
        fx_k = lam_f_x{k};
        fy_k = lam_f_y{k};
        rel_k = rsrc_k - q(k);
    
        rbm_k = [sum(fx_k); ...
                 sum(fy_k); ...
                 sum((-imag(rel_k)).*fx_k + real(rel_k).*fy_k)];
        UW((k-1)*3+1:3*k) = UW((k-1)*3+1:3*k)-rbm_k;
    end
end

%% CHECK RESIDUAL AT SURFACE
%Compute velocity at surface
B = getKmat2D(rcheck_b(1:n_bound)-q(1),0); %same for all particles
u_lhs = zeros(2*P*n_bound,1);
for k = 1:P  
    res = B*UW(3*(k-1)+1:3*k);
    u_lhs((k-1)*n_bound+1:k*n_bound) = res(1:end/2);
    u_lhs(P*n_bound+(k-1)*n_bound+1:P*n_bound+k*n_bound) = res(end/2+1:end); 
end

%Using representation
u_rhs = matvec_mob_peanut_enhanced(tau,geom_check,basis_mob);
S_0 = getRecompletionFlow(rvec_in_c,rcheck_b,q,F,T); 
u_rhs = u_rhs-S_0;  %Note! Sign here due to how we have defined the completion flow. 
                    %This is accordinng to the representation of the flow

disp('Surface residual')
diff_vec = u_rhs-u_lhs;

%max_abs = max(abs(u_rhs(1:end/2)+1i*u_rhs(end/2+1:end)));
max_abs = max(abs(S_0(1:end/2)+1i*S_0(end/2+1:end)));
res = abs(diff_vec(1:end/2)+1i*diff_vec(end/2+1:end));
abs_res = max(res); 
if max_abs > 0
    rel_vec = res/max_abs;
else
    rel_vec = res;
end
rel_res = max(rel_vec);
fprintf('Relative boundary error: %.3e\n', rel_res);
fprintf('Absolute boundary error: %.3e\n', abs_res);

 
if visualise
%     % Plot residual in x, y components separately
%     figure(9)
%     subplot(2,2,1)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(diff_vec(1:end/2))),30,log10(abs(diff_vec(1:end/2))),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('error in x velocit´y')
% 
% 
%     subplot(2,2,2)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(diff_vec(end/2+1:end))),30,log10(abs(diff_vec(end/2+1:end))),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('error in y velocity')
% 
%     % Visualise the actual velocity in the rhs and lhs
%     subplot(2,2,3)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(u_rhs(1:end/2))),30,log10(abs(u_rhs(1:end/2))),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('x velocity rhs')
% 
%     subplot(2,2,4)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(u_lhs(1:end/2))),30,log10(abs(u_lhs(1:end/2))),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('x velocity lhs')
% 
% 
%     sgtitle('Error on boundary mob peanut compression','interpreter','latex')


    rvis = [];
    aa = 0.9;
    for k = 1:P
        rvis = [rvis; aa*(rcheck_b(n_bound*(k-1)+1:k*n_bound)-q(k))+q(k)];
    end

    if strcmp(surface_error_mode,'rel')
        err_plot = rel_vec;
        err_label = 'relative';
    else
        err_plot = res;
        err_label = 'absolute';
    end
    err_plot_log = log10(max(err_plot, eps));

    figure()
    scatter3(real(rvis),imag(rvis),err_plot_log,30,err_plot_log,'filled');
    hold on
    plot(real(rcheck_b),imag(rcheck_b),'k.')
    c = colorbar;
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    c.TickLabelInterpreter = 'latex';
    ylabel(c, sprintf('%s boundary error (log10)', err_label), 'interpreter','none','FontSize',18)
    c.FontSize = 18; 
    title(sprintf('%s: Boundary %s error', solver_name, err_label), 'interpreter','none')


    %% Visualise source strengths
    figure()
    semilogy(abs(lambda_c))
    title('Source strengths mobility, peanut compression')
    
end


end


%OLD... ignore

function doPairBasisTest(Upf,Ypf,i,p2,q,U,Y,rbase_in_c,rbase_out_f,rpair_fine)
%Test pair basis with a smooth coarse density
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
        step1 = U'*[tau_p1_x;tau_p1_y]; %here I assume x and y follow each other?
        tau_mapped = Y*step1;
        Nother = singleLayer(rbase_in_c+q(i),rbase_out_f+q(p2),mu);
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
        Nother = singleLayer(rbase_in_c+q(p2),rbase_out_f+q(i),mu);
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

        Npair = singleLayer(rpair_fine,rcheck,mu);
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
        Npair = singleLayer(rpair_fine,[rbase_out_f+q(i); rbase_out_f+q(p2)],mu);
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

function test_solve_mob

close all; 
q = [0; 2.001; 2.001i]; %center coordinates

delta = 1e-3; 
P = 3; 
q = 0:2+delta:(P-1)*(2+delta);
P = 4; 
q = [0; 2+delta; 7; 9+delta];


P = 20; 
side = 2 + delta;               % neighbor center distance
R = side / (2*sin(pi/P));         % ring radius
q = R * exp(1i * (0:P-1).' * (2*pi/P));
%q(1) = 8;
%q = q+5; 
%q = q-q(1);

rng(5); 
q = grow_cluster(P,delta,2);


F = [real(q) imag(q)]; 
T = zeros(size(q)); 
rads = ones(size(q)); 
%F = [1 0; -1 0];

%F = [1 0; 0 0; 0 1; -1 0]; %forces on the particles
%T = [1; 1; 1; -1]; %torques on the particles
%rads = [1; 1; 1; 1]; 
visualise = 1; 
delta_pair = 0.5; 
%[UW1,lambdahat,it1,gmres_tol, err1] = solve_mob_2B_images(q,F,T,rads,delta_pair,visualise);

%compare to a solution with image enhancement
N_peanut = 400; 
%[UW1, lambdahat1,it1,gmres_tol, rel1, abs1] = solve_mob_2B_images(q,F,T,rads,delta_pair,visualise);
gmres_tol = 1e-8;
images = 1; 
lr = 0; 
debug = 0; 
%[UW1,lambda_mob,it1,gmres_tol,err1] = solve_mob_1B(q,F,T,rads,images, lr, visualise);
[UW2,lambdahat2,it2,gmres_tol, rel2, abs2] = solve_mob_peanut_enhanced(q,F,T,delta_pair,N_peanut,visualise,gmres_tol,debug);
% [UW3,lambdahat3,it3,gmres_tol, rel3, abs3] = solve_mob_peanut_images(q,F,T,rads,delta_pair,N_peanut,visualise);
% [UW4,lambdahat4,it4,gmres_tol, rel4, abs4] = solve_mob_2B_enhanced(q,F,T,delta_pair,visualise);
alignfigs(3);

end
