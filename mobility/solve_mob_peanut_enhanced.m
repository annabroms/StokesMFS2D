function [UW,lambda,it,gmres_tol, maxres] = solve_mob_peanut_enhanced(q,F,T,delta_pair,N_peanut,visualise)
%SOLVE_MOB_PEANUT_ENHANCED Solve a 2D Stokes mobility problem with peanut-compressed pair corrections.
%
% Syntax:
%   [UW,lambda,it,gmres_tol,maxres] = solve_mob_peanut_enhanced( ...
%       q,F,T,rads,delta_pair,N_peanut,visualise)
%
% Inputs:
%   q          - Complex particle centers (P x 1).
%   F          - Net forces per particle (P x 2), columns are x/y.
%   T          - Net torques per particle (P x 1).
%   rads       - Particle radii (P x 1).
%   delta_pair - Near-pair threshold; pairs closer than this use local fine correction.
%   N_peanut   - Number of peanut check points used in pair compression.
%   visualise  - Plot diagnostics if true.
%
% Outputs:
%   UW         - Rigid-body motion vector [Ux1; Uy1; W1; ...; UxP; UyP; WP].
%   lambda     - Stacked coarse source strengths [lambda_x; lambda_y].
%   it         - GMRES iterations.
%   gmres_tol  - GMRES tolerance used.
%   maxres     - Maximum relative boundary residual on an independent check grid.
%
% Notes:
%   - Call with no inputs to run the built-in test.
%   - Pair-source contributions in UW include both body-fine and image-fine nodes.
%
% Anna Broms, Feb 23, 2026

if nargin==0, test_solve_mob; 
    return; end


%% SET PARAMS
%GMRES params
maxit = 800; 
gmres_tol = 1e-6; %not enough given the residual we seek
gmres_tol = 1e-10; 

% Grid params
P = length(q); 
%Set coarse and fine grid. 

%Play with N_c, N_f, a (a_f). 
N_c = 60;  %100 better here? 
%N_c = 100; 
%N_c = 200; 
N_f = 60; 
%N_f = N_c; %debug

%N_c = 250; 
a_c = 1.2;
%a_c = 2.5;
%a_c = 5;
 

a_f = 1.2; %upsampling factor for the fine grid

tol_c = 1e-12; %I think this works reasonably
%tol_c = 1e-10; %Curve moves closer to the surface -> smaller coeff  
%tol_c = 1e-16; %Curve moves further from surface -> larger coeff. 

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rp_c = max([1-sep_c,0.01]); %radius of proxy surface for coarse grid
Rp_f = max([1-sep_f,0.01]);  % and fine grid
%Rp_c = 0.85;
%Rp_f = Rp_c; %debug

%accumulation point, given Rg and delta. Closed formula from fixed point of reflection formula
accstop = (1-Rp_c)^2/Rp_c;  

if nargin < 6
    delta_pair = accstop; %We want to use the pair correction for all gaps smaller than delta_pair. (or accstop).
end

opt = get2Dparams(); 
opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c; 
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 
opt.N_peanut = N_peanut;
opt.precomp = 1; %faster if evaluation of one body basis on fine grid is compted only once. 
% %Less storage required.
opt.pc = 1; %prepare grid to do pair corrections
opt.delta_pair = delta_pair; 

opt.cmap = 1; 
opt.Nclust = 200; %points on ellipse segments, for now
 

%% CREATE GRID
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


%get evaluation of lambda0
[urhs,~] = getRecompletionFlow(rvec_in_c,rout,q,F,T); 

Kf = getKmat2D(rbase_in_f,0);
Lf = Kf*((Kf'*Kf)\Kf'); %This is x y

%Get one-body pseduo inverse blocks -- enough to do this for single body.
rimage_in = []; 
[U,Y,Lc] = getSelfPseudoMobilityStokes(1,q,rbase_in_c,rbase_out_c,rimage_in,[0,ceil(a_c*N_c)]);

%Get pair basis
[UB_all,YB_all,UC_all,YC_all,Cmap,~] = getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,1,Lc{1});
                                    
%build the mapping that maps coarse boundary data to coarse proxy data
%A = buildCoarseMap(UB_all,YB_all,UC_all,YC_all,U,Y);


%Visualisation for presentation / paper
% if visualise
%     figure(2) 
%     hold on
%     
%     tout_f = linspace(0,2*pi,ceil(a_f*N_f)+1);
%     tout_f = tout_f(1:end-1)';
%     rbase_out_f = cos(tout_f)+1i*sin(tout_f);
%     rin_f = [];
%     rout_plot = [];
%     for k = 1:P
%         rin_f = [rin_f; q(k)+rbase_in_f]; %just for visuals
%         rout_plot = [rout_plot; q(k)+rbase_out_c]; %just for visuals
%     end
%     fill(real(rout_plot(1:end/2)),imag(rout_plot(1:end/2)),[0 0 0],'FaceAlpha',0.1,'EdgeColor','none')
%     hold on
%     fill(real(rout_plot(end/2+1:end)),imag(rout_plot(end/2+1:end)),[0 0 0],'FaceAlpha',0.1,'EdgeColor','none')
% 
%     red = [231,51,57]/256;
%     yellow = [251,186,0]/256;
%     blue = [25,105,188]/256;
%     newred = 0.7*red+0.3*yellow;
%     plot(real(rvec_in_c),imag(rvec_in_c),'k.');
%     hold on
%     %plot(real(rin_f),imag(rin_f),'.','Color',newred);
%     
%     
% %     rf1 = refine{1,2};
% %     rf2 = refine{2,1};
% %     plot(real(rout_f),imag(rout_f),'.','Color',blue)
% %     plot(cos(rf1),sin(rf1),'.','Color',blue)
% %     plot(cos(rf2)+real(q(2)),sin(rf2)+imag(q(2)),'.','Color',blue)
%     
%     plot(real(rout_plot),imag(rout_plot),'.','Color',blue)
%     rimage1 = rimage_vec{2,1};
%     rimage2 = rimage_vec{1,2};
%     rimage1 = rimage1(end);
%     rimage2 = rimage2(1); 
%     plot(real(rimage1),imag(rimage1),'.','Color',newred,'MarkerSize',18)
%     plot(real(rimage2),imag(rimage2),'.','Color',newred,'MarkerSize',18);
%     axis equal
% 
%     delta = abs(q(1)-q(2))-2;
%     alpha = acos((2+delta)/4);
%     x = 1+delta/2;
%     y = 2*sin(alpha);
% 
% %     third = x+cos(tout_c_all)+1i*y+1i*sin(tout_c_all);
% %     plot(real(third),imag(third),'k--');
% % 
% %     third2 = x+cos(tout_c_all)-1i*y+1i*sin(tout_c_all);
% %     plot(real(third2),imag(third2),'k--');
% 
% end




%Visualise 1-body and pair-basis
%warning('Deactivate opt.precomp');
%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,[],[],N_c, N_f,a_c,a_f,rads)


Lc_pair = getILpair(Lc{1});
Lf_pair = getILpair(Lf);


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
basis_mob.Lc_pair = Lc_pair;
basis_mob.Lf_pair = Lf_pair;

geom_check = geom_solve;
geom_check.rcheck = rcheck_b;

%% Solve system

% Build the matrix to check it out
debug = 1;
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
    skeel(CC)

    [V,D] = eig(CC);
    D = diag(D); 
    figure()
    plot(real(D),imag(D),'+')
    xlabel('Re \lambda')
    ylabel('Im \lambda')
end

debug = 0;                                                                                                                                           
[tau,it,resvec,real_res] = helsing_gmres(@(x) matvec_mob_peanut_enhanced(x,geom_solve,basis_mob),...
    urhs,2*size(rout,1),maxit,gmres_tol,1,rout);

figure()
semilogy(resvec); 
title('GMRES convergence with peanut compression, mobility', 'interpreter','latex')


if debug
    %check residual
    restot = (matvec_mob_pairprecond_peanut(tau,rbase_in_c,rbase_in_f,rvec_in_c,refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair)-urhs)./urhs;
    figure()
    semilogy(abs(restot))
    title('Rel res at colloc points for mob peanut')
end
% hold on
% semilogy(resvec2); 

%% COMPUTE Rigid body motion
%And evaluate residual in new points rcheck_b

[tau_stokes_x, tau_stokes_nonpx,~, tau_beta_x,tau_stokes_y, ...
    tau_stokes_nonpy,~,tau_beta_y,~,rimage_k] = ...
    transform_mob_peanut_stokes(tau,geom_check,basis_mob);
lambda = [tau_stokes_x; tau_stokes_y];

%%% Get rigid body motion. 

%First due to all coarse sources
Kc = getKmat2D(rbase_in_c,0);
UW= zeros(3*P,1); 
for k= 1:P
    UW((k-1)*3+1:3*k) = -Kc'*[tau_stokes_nonpx((k-1)*N_c+1:k*N_c); tau_stokes_nonpy((k-1)*N_c+1:k*N_c)];
end

% Then, due to all pair sources (fine-body + fine-image).
has_neigh = sort(unique(pairs(:)));
for i = 1:length(has_neigh)
    k = has_neigh(i);
    rsrc_k = [rbase_in_f+q(k); rimage_k{k}];
    fx_k = tau_beta_x{k};
    fy_k = tau_beta_y{k};
    rel_k = rsrc_k - q(k);

    rbm_k = [sum(fx_k); ...
             sum(fy_k); ...
             sum((-imag(rel_k)).*fx_k + real(rel_k).*fy_k)];
    UW((k-1)*3+1:3*k) = UW((k-1)*3+1:3*k)-rbm_k;
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

%debug                    
% figure()
% subplot(1,2,1)
% plot(u_rhs);
% hold on
% plot(u_lhs);
% subplot(1,2,2)
% semilogy(abs(u_rhs-u_lhs));


disp('Surface residual')
diff_vec = u_rhs-u_lhs;
max_abs = max(abs(u_rhs(1:end/2)+1i*u_rhs(end/2+1:end)));
err_vec =  abs(diff_vec(1:end/2)+1i*diff_vec(end/2+1:end))/max_abs;
maxres = max(err_vec)


 
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
%     title('error in x velocity')
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

    figure()
    scatter3(real(rvis),imag(rvis),log10(abs(err_vec)),30,log10(abs(err_vec)),'filled');
    hold on
    plot(real(rcheck_b),imag(rcheck_b),'k.')
    c = colorbar;
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    c.TickLabelInterpreter = 'latex';
    ylabel(c,'Relative residual $\log_{10}$','interpreter','latex','FontSize',18)
    c.FontSize = 18; 
    title('Mobility peanut')


    %% Visualise source strengths
%     figure()
%     semilogy(abs(lambda))
%     title('Source strengths mobility, peanut compression')
    
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

delta = 1e-2; 
P = 3; 
q = 0:2+delta:(P-1)*(2+delta);
P = 4; 
q = [0; 2+delta; 5; 7+delta];


P = 2; 
side = 2 + delta;               % neighbor center distance
R = side / (2*sin(pi/P));         % ring radius
q = R * exp(1i * (0:P-1).' * (2*pi/P));
%q(1) = 8;
F = [real(q) imag(q)]; 
T = zeros(size(q)); 
rads = ones(size(q)); 

%F = [1 0; 0 0; 0 1; -1 0]; %forces on the particles
%T = [1; 1; 1; -1]; %torques on the particles
%rads = [1; 1; 1; 1]; 
visualise = 1; 
delta_pair = 0.2; 
%[UW1,lambdahat,it1,gmres_tol, err1] = solve_mob_precond_images(q,F,T,rads,delta_pair,visualise);

%compare to a solution with image enhancement
N_peanut = 400; 
[UW2,lambdahat,it2,gmres_tol, err2] = solve_mob_peanut_enhanced(q,F,T,delta_pair,N_peanut,visualise);
[UW3,lambdahat,it3,gmres_tol, err3] = solve_mob_precond_peanut(q,F,T,rads,delta_pair,N_peanut,visualise);

str = sprintf('Relative residual with 2-body precond: %1.2e vs with peanut compression: %1.2e\n Converging in %u resp % u iterations',err1,err2,it1,it2);
disp(str)

alignfigs(3);

end
