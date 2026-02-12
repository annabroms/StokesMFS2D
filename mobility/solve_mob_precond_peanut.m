function [UW,lambda,it,gmres_tol, maxres] = solve_mob_precond_peanut(q,F,T,rads,delta_pair,N_peanut,lr,visualise)
%SOLVE_MOB_PRECOND_PEANUT Solves a 2D Stokes mobility problem with circular
%particles using a 2-body preconditioned recompleted MFS formulation. A
%fine grid enhanced with approximate image points is used locally for every
%close pair of particles to resolve challenging interactions. This fine grid is
%compressed so that only coarse sources and collocation points are needed to 
%solve the problem iteratively, effectively preconditioning the system.
%
% Syntax:
%   [UW, lambdahat, it, gmres_tol, maxres] = solve_mob_precond_peanut(q, F, T, rads, delta_pair, N_peanut, visualise)
%
% Inputs:
%   q          - Vector of length P, complex-valued center coordinates for the particles
%   F          - Px2 matrix of net force vectors (columns: x and y components)
%   T          - Px1 column vector of torques acting on the particles
%   rads       - Px1 vector of particle radii
%   delta_pair - Scalar threshold used to determine which particle pairs are considered close. For such pairs, a fine BVP is solved locally (a pair correction is built).
%   N_peanut   - Number of points on the peanut separation surface between
%                every close pair of particles. The peanut boundary is used
%                to map fine sources to effective coarse sources, giving the
%                same flow field exterior to the close pair of particles.
%   lr         - Flag for long range preconditioning 
%   visualise  - Logical flag: plot the configuration and solution details

%
% Outputs:
%   UW         - 3P×1 vector of computed rigid-body motion (RBM) velocities
%   lambda     - Solution vector of source strengths
%   it         - Number of GMRES iterations required
%   gmres_tol  - Set GMRES tolerance
%   maxres     - Maximum relative residual in a test (non-collocation) set of boundary nodes
%
% Description:
%   The FMM is used for Stokeslet evaluation. No other source types are
%   needed thanks to compression of the local fine grid.
%
% Notes:
%   - Intended to match solve_precond_peanut, but applied to a mobility problem rather than resistance.
%   - Aims to test an MFS generalisation of the idea presented by Cheng-Greengard (1998)
%
% See also:
%   solve_2D_mob              - 1-body preconditioned mobility solver
%   solve_mob_precond_images - 2-body preconditioner without peanut compression
%   solve_precond_peanut   - 2-body preconditioned resistance solver
%   with peanut compression
%
% To test: Call without arguments.
%
% Anna Broms, April 4, 2025

if nargin==0, test_solve_mob; 
    return; end


%% SET PARAMS
%GMRES params
maxit = 800; 
gmres_tol = 1e-6; %not enough given the residual we seek
gmres_tol = 1e-8; 

% Grid params
P = length(q); 
%Set coarse and fine grid. 

%Play with N_c, N_f, a (a_f). 
N_c = 150;  %100 better here? 
N_c = 80; 
%N_c = 200; 
N_f = 150; 
%N_f = N_c; %debug

%N_c = 250; 
a_c = 1.2;
%a_c = 2.5;
%a_c = 5;
 

a_f = 1.2; %upsampling factor for the fine grid

tol_c = 1e-12; %I think this works reasonably
%tol_c = 1e-10; %Curve moves closer to the surface -> smaller coeff  
%tol_c = 1e-16; %Curve moves further from surface -> larger coeff. 

s = [0 0 1 1 0 0 0];  %set type of singularities at image points [S R T D]
%s = [1 0 1 1]; %Other singularities? Currently not supported! But code can
%be changed!

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
opt.s = s;
opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.P = P; 
opt.a_c = a_c; 
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 
opt.lr = lr;
opt.rads = rads; 
opt.N_peanut = N_peanut;
opt.s = s; 
opt.precomp = 1; %faster if evaluation of one body basis on fine grid is compted only once. 
% %Less storage required.
opt.pc = 1; %prepare grid to do pair corrections
opt.delta_pair = delta_pair; 
opt.cmap = 1; 
 

%% CREATE GRID
%Outer basic grid
tout_c_all = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c_all(1:end-1)';
rbase_out_c = cos(tout_c)+1i*sin(tout_c);
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rp_c*cos(tin)+Rp_c*1i*sin(tin); 
rvec_in_c = [];
for k = 1:P
    rvec_in_c = [rvec_in_c; q(k)+rbase_in_c];
end

%Construct image grid
%[rout, weights, rin, rimage, nimage, pair_points, pairs, rimage_pairs, refine, rin_base] 
[rout,~,~,~,~,~,pairs,rimage_vec,refine,rbase_in_f] = get2DImageGrid(q,rads,opt);


%get evaluation of lambda0
[urhs,~] = getRecompletionFlow(rvec_in_c,rout,q,F,T); 


Kf = getKmat2D(rbase_in_f,0);
Lf = Kf*((Kf'*Kf)\Kf'); %This is x y


%Get one-body pseduo inverse blocks -- enough to do this for single body.
[U,Y,Lc,Lr] = getSelfPseudoMobility(1,q,rbase_in_c,rbase_out_c,[],[],[0,ceil(a_c*N_c)],s,1);


%Get pair basis
[UB_all,YB_all,UC_all,YC_all,Cmap,~,nimage] = getPairBasis(q,N_f,a_f,rads,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,1,Lc{1},Lf,Kf);

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
    rcheck_b = [rcheck_b; q(k)+rads(k)*(cos(t)+1i*sin(t))];
end



%% Experiment with left preconditioner based on deflation
debug = 1;
if debug && lr
    x = zeros(2*length(rout),1);
    tic
    for k = 1:2*length(rout)
        k
        x(:) = 0; 
        x(k) = 1; 

        uu = matvec_mob_pairprecond_peanut(x,rbase_in_c,rbase_in_f,rvec_in_c,...
             refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all,YC_all,Cmap,Lc_pair,Lf_pair);
      
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

    num_eigs = 3; 
    [ss,I] = maxk(abs(D),num_eigs);
    Vk = V(:,I).*ss';

    [Uc,S,V] = svd(CC); 
    Vsmall = V(:,end); 
    

    Mc = round(opt.a_c*opt.N_c);
    t = linspace(0,2*pi,Mc+1);
    t = t(1:end-1)';  
    for i = 1
        for k = 1:P  
            Vpx = Vsmall((k-1)*Mc+1:k*Mc,i);
            Vpy = Vsmall((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
    
            figure(6)
           % subplot(1,3,i)
            quiver(real(q(k))+cos(t),imag(q(k))+sin(t),Vpx,Vpy);
            hold on
            colorbar
            view(0,90)
        end
    end





end

if lr
    rin_c = []; 
    for k = 1:P
        rin_c = [rin_c; rbase_in_c+q(k)];
    end
                     %get_long_range_precond_peanut_mob(q,rin,rout,L,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair,opt)

    [Sinv,Zi,Yi,db] = get_long_range_precond_peanut_mob(q,rin_c,rout,Lc{1},Lr,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,U,Y,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair,opt,Vsmall);

    %[Sinv,Zi,Yi,db] = get_long_range_precond_mu(q,rin_c,rout,opt);
    opt.db = db;
    mu_coarse = getCoarseMu(urhs,Sinv,Zi,Yi,db,P,opt.N_c,opt.a_c,Vsmall);
   % mu_coarse = Vsmall*Sinv*Vsmall'*urhs;
    %mu_coarse = Vsmall*Vsmall'*urhs; %Better to do the projection like
    %this? 

    %Try to evaluate: 
    res1 = matvec_mob_pairprecond_peanut(mu_coarse,rbase_in_c,rbase_in_f,rvec_in_c,refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair);
    
    %pair_points = [zeros(P,1) opt.a_c*opt.N_c*ones(P,1)];
    %res2 = matvec_2D_Stokes(mu_coarse,rvec_in_c,rout,[],[],q,UU,YY,pair_points,s);
end

%% Solve system

% Build the matrix to check it out
debug = 0;
if debug
    x = zeros(2*length(rout),1);
    tic
    for k = 1:2*length(rout)
        k
        x(:) = 0; 
        x(k) = 1; 
        if lr 
            uu = lr_matvec_mob_peanut(x,rin_c,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair,Sinv,Zi,Yi);
        else
            uu = matvec_mob_pairprecond_peanut(x,rbase_in_c,rbase_in_f,rvec_in_c,...
                refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all,YC_all,Cmap,Lc_pair,Lf_pair);
        end
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

if lr
    Pf = applyPmat_peanut_mob(urhs,rin_c,rout,Lc{1},Sinv,q,Zi,Yi,rbase_in_c,...
        rbase_in_f,refine,rimage_vec,nimage,opt,U,Y,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair,Vsmall);   
   [tau,it,resvec,real_res] = helsing_gmres(@(x) lr_matvec_mob_peanut(x,rin_c,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair,Sinv,Zi,Yi,Vsmall),Pf,2*size(rout,1),maxit,gmres_tol,1,rout);   

else

    debug = 0;                                                                                                                                           
    [tau,it,resvec,real_res] = helsing_gmres(@(x) matvec_mob_pairprecond_peanut(x,rbase_in_c,rbase_in_f,rvec_in_c,refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair),urhs,2*size(rout,1),maxit,gmres_tol,1,rout);
end
figure()
semilogy(resvec); 
title('GMRES convergence with peanut compression, mobility', 'interpreter','latex')


if visualise
    %check residual
    if lr
        restot = (lr_matvec_mob_peanut(tau,rin_c,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair,Sinv,Zi,Yi,Vsmall)-Pf)./Pf;
    else
        restot = (matvec_mob_pairprecond_peanut(tau,rbase_in_c,rbase_in_f,rvec_in_c,refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair)-urhs)./urhs;
    end
    figure()
    semilogy(abs(restot))
    title('Rel res at colloc points for mob peanut')
end
% hold on
% semilogy(resvec2); 

%% COMPUTE Rigid body motion
%And evaluate residual in new points rcheck_b

if lr
    mu_mapped = applyQmat_peanut_mob(tau,rvec_in_c,rout,Sinv,Zi,Yi,opt,...
            rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair);
    
    tau = mu_mapped+mu_coarse; 
end

[tau_stokes_x, tau_stokes_nonpx,~, tau_beta_x,tau_stokes_y,tau_stokes_nonpy,~,tau_beta_y,~] = transform_mob_peanut(tau,...
    rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,rcheck_b,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair);
lambda = [tau_stokes_x; tau_stokes_y];



%%% Get rigid body motion. 

%First due to all coarse sources
Kc = getKmat2D(rbase_in_c,0);
UW= zeros(3*P,1); 
for k= 1:P
    UW((k-1)*3+1:3*k) = -Kc'*[tau_stokes_nonpx((k-1)*N_c+1:k*N_c); tau_stokes_nonpy((k-1)*N_c+1:k*N_c)];
end

% %Then, due to all fine sources.
has_neigh = sort(unique(pairs(:)));
for i = 1:length(has_neigh)
    k = has_neigh(i); 
    UW((k-1)*3+1:3*k) = UW((k-1)*3+1:3*k)-Kf'*[tau_beta_x((k-1)*N_f+1:k*N_f); 
        tau_beta_y((k-1)*N_f+1:k*N_f)];
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
u_rhs = matvec_mob_pairprecond_peanut(tau,rbase_in_c,rbase_in_f,rvec_in_c,refine,rimage_vec,nimage,opt,rout,rcheck_b,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair);
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
P = 4; 
%alternative configuration
delta = 0.001; 
x = 1+delta/2;
y = sqrt((2+delta)^2-(1+delta/2)^2);
q = [0; 2+delta; x+1i*y];
%q = hexagonal_lattice(delta, 2);
% q = q(1:3); 
%F = [1 0; 0 0; 0 1]; %forces on the particles
%T = [1; 1; 1]; %torques on the particles
P = length(q);
%q = [0:(2+delta):(P-1)*(2+delta)]';
F = rand(P,2); T = rand(P,1); 
rads = ones(P,1); 
visualise = 1; 
delta_pair = 0.2;

%delta_pair = 2; 
%[UW1,lambdahat,it1,gmres_tol, err1] = solve_mob_precond_images(q,F,T,rads,delta_pair,visualise);

%compare to a solution with image enhancement
N_peanut = 400;
lr = 20; 
%[UW2,lambdahat,it2,gmres_tol, err2] = solve_mob_precond_peanut(q,F,T,rads,delta_pair,N_peanut,lr,visualise);
[UW3,lambdahat,it3,gmres_tol, err3] = solve_mob_precond_peanut(q,F,T,rads,delta_pair,N_peanut,0,visualise);
%str = sprintf('Relative residual with 2-body precond: %1.2e vs with peanut compression: %1.2e\n Converging in %u resp % u iterations',err1,err2,it1,it2);
%disp(str)

alignfigs(3);

end



