function [UW, lambda, it, gmres_tol, maxres] = solve_mob_precond_enhanced(q, F, T, rads, delta_pair, visualise)
%SOLVE_MOB_PRECOND_ENHANCED Solves a 2D Stokes mobility problem with circular
%particles using a 2-body preconditioned recompleted MFS formulation. To resolve
% challenging close interactions, a fine 2-body BVP is solved for fine
% source strengths on each close pair of particles. The fine sources,
% correct the representation obtained from a coarse grid, effectively preconditioning the system.
%
% Syntax:
%   [UW, lambda, it, gmres_tol, maxres] = solve_mob_precond_enhanced(q, F, T, rads, delta_pair, visualise)
%
% Inputs:
%   q          - Vector of length P, complex-valued center coordinates for the particles
%   F          - Px2 matrix of net force vectors (columns: x and y components)
%   T          - Px1 column vector of torques acting on the particles
%   rads       - Px1 vector of particle radii
%   delta_pair - Scalar threshold used to determine which particle pairs are considered close. For such pairs, a fine BVP is solved locally (a pair correction is built).
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
%   This function applies a 2-body preconditioner (using pair corrections via local fine BVPs)
%   along with image systems to improve solution accuracy for near-contact particle configurations.
%   The FMM is used for Stokeslet evaluation, while other source types use direct summation.
%
% Notes:
%   - Intended to match solve_2D_precond_images, but applied to a mobility problem rather than resistance.
%
% See also:
%   solve_2D_mob              - 1-body preconditioned mobility solver
%   solve_mob_precond_peanut - 2-body preconditioner with peanut
%                              compression (in this version of the algorithm, only coarse sources
%                              are needed at the solve stage).
%   solve_res_precond_enhanced   - 2-body preconditioned resistance solver
%
% To test: Call without arguments.
%
% Anna Broms, Feb12, 2026

if nargin==0, test_solve_mob; 
    return; end


%% SET PARAMS
%GMRES params
maxit = 800; 
gmres_tol = 1e-6; 
gmres_tol = 1e-10; 

%Grid params
P = length(q); 

opt = get2Dparams();
opt.P = P; 

%Set coarse and fine grid. 
%Play with N_c, N_f, a (a_f). 
N_c = 60;  %100 better here? 
%N_c = 60; 
%N_c = 100;  
N_f = 60; 


a_c = 1.2;
%a_c = 1; 

 
a_f = 1.2; %upsampling factor for the fine grid
%a_f = 5;  

tol_c = 1e-12; %I think this works reasonably
%tol_c = 1e-8; %Curve moves closer to the surface -> smaller coeff and 
% smoother coarse 1-body basis to evaluate on neighbour
%tol_c = 1e-16; %Curve moves further from surface -> larger coeff. 

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rp_c = max([1-sep_c,0.01]); %radius of proxy surface for coarse grid
Rp_f = max([1-sep_f,0.01]);  % and fine grid
%Rp_f = Rp_c; %debug

%accumulation point, given Rp and delta. Closed formula from fixed point of reflection formula
accstop = (1-Rp_c)^2/Rp_c;  

if nargin < 5
    delta_pair = accstop; %We want to use the pair correction for all gaps smaller than delta_pair. (or accstop).
end


opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c; 
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 
opt.N_peanut = 0;
opt.precomp = 1; %faster if evaluation of one body basis on fine grid is compted only once. 
% %Less storage required.
opt.pc = 1; %prepare grid to do pair corrections
opt.delta_pair = delta_pair; 
opt.Nclust = 200;
opt.use_cached_pair_transform = false; % set true to use getMobPairTransformationStokesCached


%% CREATE GRID
%Outer basic grid
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';

rbase_out_c = rads(1)*cos(tout_c)+1i*rads(1)*sin(tout_c);
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rp_c*cos(tin)+Rp_c*1i*sin(tin);

tin_f = linspace(0,2*pi,N_f+1);
tin_f = tin_f(1:end-1)';
rbase_in_f =  Rp_f*cos(tin_f)+Rp_f*1i*sin(tin_f);
 
[~, ~, ~, rimage_vec, refine,pairs] = getEnhancedGrid(q, opt);

rin = [];
rout = [];

for k = 1:P
    rin = [rin; q(k)+rbase_in_c];
    rout = [rout; rbase_out_c+q(k)];
end

% % To illustrate for presentation / paper
% for k = 1:P
%     rin_f = [rin_f; q(k)+rbase_in_f]; %just for visuals
%     rout_f = [rout_f; q(k)+rbase_out_f]; %just for visuals
% end
% if visualise
%     %Visualise the compression between coarse and fine grid
%     figure()
%     fill(real(rout_f(1:end/2)),imag(rout_f(1:end/2)),[0 0 0],'FaceAlpha',0.1,'EdgeColor','none')
%     hold on
%     fill(real(rout_f(end/2+1:end)),imag(rout_f(end/2+1:end)),[0 0 0],'FaceAlpha',0.1,'EdgeColor','none')
%     %set colors for presentation... 
%     red = [231,51,57]/256;
%     yellow = [251,186,0]/256;
%     blue = [25,105,189]/256;
%     newred = 0.7*red+0.3*yellow;
%     plot(real(rin),imag(rin),'k.');
%     hold on
%     %plot(real(rin_f),imag(rin_f),'.','Color',newred);
%     
%     
%     rf1 = refine{1,2};
%     rf2 = refine{2,1};
%     %plot(real(rout_f),imag(rout_f),'.','Color',blue)
%     %plot(cos(rf1),sin(rf1),'.','Color',blue)
%     %plot(cos(rf2)+real(q(2)),sin(rf2)+imag(q(2)),'.','Color',blue)
%     
%     plot(real(rout),imag(rout),'b.')
%     rimage1 = rimage_vec{2,1};
%     rimage2 = rimage_vec{1,2};
%     %plot(real(rimage1),imag(rimage1),'.','Color',newred);
%     %plot(real(rimage2),imag(rimage2),'.','Color',newred);
%     axis equal
% end

%% PREPARE PRECONDITIONING AND RHS

%get evaluation of lambda0, the completion sources, computed from known net forces and torques
[urhs,lambda_x,lambda_y] = getRecompletionFlow(rin,rout,q,F,T); 

% The pseudo inverse for a single particle is supposed to use the 
% coarse proxy grid + basic grid of collocation points only. We assume
% everyone has the same discretisation, so it's enough to do this once.
rimage_in = []; 
[U,Y,Lc] = getSelfPseudoMobilityStokes(1,q,rbase_in_c,rbase_out_c,rimage_in,[0,ceil(a_c*N_c)]);

%Kf = getKmat2D(rbase_in_f,0); % Kf' maps force density to net force and torque
%Lf = Kf*((Kf'*Kf)\Kf'); %Projects onto the range of the constraint matrix Kf'

%Get pair basis
%[Upf,Ypf,~,~,~,~,nimage] = getPairBasis(q,N_f,a_f,rads,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,1,Lc{1},Lf,Kf);
%Get pair basis
[Upf,Ypf,~,~,~,~] = getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,1,Lc{1});

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.rimage_vec = rimage_vec;
geom.opt = opt;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Lc = Lc{1};
basis.Upf = Upf;
basis.Ypf = Ypf;


% Now, check pair basis up to the boundary. Is it nice and smooth?
%warning('Deactivate opt.precomp');
%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,Lc{1},Lf,N_c,N_f,a_c,a_f,rads)

%checkOneBasis(rbase_in_c,U,Y,Lc{1},N_c,a_c)
%% Construct check boundaries
% Create new grid points, in which the accuracy of the solution is
% to be evaluated. 
rcheck_b = [];
n_bound = 803;
t = linspace(0,2*pi,n_bound)';
for k = 1:P
    rcheck_b = [rcheck_b; q(k)+rads(k)*(cos(t)+1i*sin(t))];
end

%% SOLVE SYSTEM
% Build the matrix to check it out

debug = 1;
if debug
    x = zeros(size(urhs));
    for k = 1:size(urhs,1)
        k
        x(:) = 0; 
        x(k) = 1; 
        uu = matvec_mob_pairprecond_enhanced(x,geom,basis);
        CC(:,k) = uu;
    end
    figure()
    imagesc(log10(abs(CC)))
    colorbar
    skeel(CC)

    get_nullspace = 0;
    get_eigs = 1; 

    if get_eigs

        [V,D] = eig(CC);
        D = diag(D); 
        figure()
        plot(real(D),imag(D),'+')
        xlabel('Re \lambda')
        ylabel('Im \lambda')
            
        [s,I] = mink(abs(D),3);

        if get_nullspace
            Vsmall = V(:,I).*s';
        
            [UU,S,VV] = svd(CC);
            SS = diag(S);
            figure()
            semilogy(SS)
            semilogy(SS,'+')
            [s2,I2] = mink(SS,6);
        
            UUs = UU(:,I2);
            VVs = VV(:,I2); 
        
            Mc = round(a_c*N_c);
            t = linspace(0,2*pi,Mc+1);
            t = t(1:end-1)'; 
        
            %Visualise left and right singular vectors
            for i = 1:6
                for k = 1:P
                    Vpx = VVs((k-1)*Mc+1:k*Mc,i);
                    Vpy = VVs((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
                    Vp = abs(Vpx+1i*Vpy);
                    figure(33)
                    subplot(2,3,i)
                    scatter3(real(q(k))+cos(t),imag(q(k))+sin(t),Vp,40,Vp,'filled');
                    hold on
                    colorbar
                    view(0,90)
                    sgtitle('Right sing vec','interpreter','latex')
                    title(i)
                    axis off
                    axis square
        
                    figure(34)
                    subplot(2,3,i)
                    quiver(real(q(k))+cos(t),imag(q(k))+sin(t),Vpx,Vpy)
                    hold on
                    q2 = [q; q(1)];
                    plot(real(q2),imag(q2),'k-')
                    axis off
                    sgtitle('Right sing vec','interpreter','latex')
                    title(i)
                    axis square
        
        
                    Vpx = UUs((k-1)*Mc+1:k*Mc,i);
                    Vpy = UUs((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
                    Vp = abs(Vpx+1i*Vpy);
                    figure(35)
                    subplot(2,3,i)
                    scatter3(real(q(k))+cos(t),imag(q(k))+sin(t),Vp,40,Vp,'filled');
                    hold on
                    colorbar
                    view(0,90)
                    sgtitle('Left sing vec','interpreter','latex')
                    title(i)
                    axis off
                    axis square
        
                    figure(36)
                    subplot(2,3,i)
                    quiver(real(q(k))+cos(t),imag(q(k))+sin(t),Vpx,Vpy)
                    hold on
                    q2 = [q; q(1)];
                    plot(real(q2),imag(q2),'k-')
                    axis off
                    sgtitle('Left sing vec','interpreter','latex')
                    title(i)
                    axis square
                end
            end
                    
            %Visualise eigenvectors
        
            for i = 1:3
             
                for k = 1:P  
                    Vpx = Vsmall((k-1)*Mc+1:k*Mc,i);
                    Vpy = Vsmall((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
                    Vp = abs(Vpx+1i*Vpy);
                    figure(31)
                    subplot(1,3,i)
                    scatter3(real(q(k))+cos(t),imag(q(k))+sin(t),Vp,40,Vp,'filled');
                    hold on
                    colorbar
                    view(0,90)
                    sgtitle('Abs of eigvec','interpreter','latex')
                    axis off
                    axis square
        
        
                    figure(32)
                    subplot(1,3,i)
                    quiver(real(q(k))+cos(t),imag(q(k))+sin(t),Vpx,Vpy)
                    hold on
                    q2 = [q; q(1)];
                    plot(real(q2),imag(q2),'k-')
                    sgtitle('Eigvec','interpreter','latex')
                    axis off
                    axis square
                end
        
    
            end
        end

    end
    

end

[tau,it,resvec,real_res] = helsing_gmres(@(x) matvec_mob_pairprecond_enhanced(x,geom,basis),urhs,2*length(rout),maxit,gmres_tol,1,rout);
debug = 1; 

%Modify to build with krylov preconditioning
%[tau, e2, precond] = precond_gmres(@(x) ...)
%it = length(e2); 

if debug
      figure()
%     semilogy(e2);
      semilogy(resvec);
      title('GMRES convergence mobility, pair corrections', 'Interpreter','latex')
      u2 = matvec_mob_pairprecond_enhanced(tau,geom,basis);
end 

%% POSTPROCESS
%[rvec_in,rimage_in,nimage_in,coarse_ind,tau_stokes_x,tau_stokes_y, ...
 %   tau_stokes_nonpx, tau_stokes_nonpy,tau_stress_x,tau_stress_y,tau_stress_all_x,tau_stress_all_y,tau_stress_all_px,tau_stress_all_py,...
 %   tau_pot_x,tau_pot_y,tau_pot_all_x,tau_pot_all_y,tau_pot_all_px,tau_pot_all_py,u_corr]

%warning('work to be done for postprocessing here!')

if isfield(opt,'use_cached_pair_transform') && opt.use_cached_pair_transform
    [~,~,tau_stokes_x,tau_stokes_y, ...
        tau_stokes_nonpx, tau_stokes_nonpy,tau_stokes_e_nonpx, tau_stokes_e_nonpy, rimage_k, u_corr] = ...
        getMobPairTransformationStokesCached(tau,geom,basis);
else
    [~,~,tau_stokes_x,tau_stokes_y, ...
        tau_stokes_nonpx, tau_stokes_nonpy,tau_stokes_e_nonpx, tau_stokes_e_nonpy, rimage_k, u_corr] = ...
        getMobPairTransformationStokes(tau,geom,basis);
end

lambda = [tau_stokes_x; tau_stokes_y]; %This is the fine density
%tau_proxy = [tau_stokes_nonpx; tau_stokes_nonpy];
%And evaluate in new points rcheck_dom and rcheck_b

%The total source strengths are obtained by adding on lambda_x and lambda_y
%here. 

%lambdahat = [tau_proxy; tau_image]; %not just for other things than visualisation

%% Get rigid body motion. 

%First due to all coarse sources
Kc = getKmat2D(rbase_in_c,0);
UW= zeros(3*P,1); 
for k= 1:P
    UW((k-1)*3+1:3*k) = -Kc'*[tau_stokes_nonpx((k-1)*N_c+1:k*N_c); tau_stokes_nonpy((k-1)*N_c+1:k*N_c)];
end

%Then, due to all fine sources.
has_neigh = sort(unique(pairs(:)));
Kf = getKmat2D(rbase_in_f,0); 
for i = 1:length(has_neigh)
    k = has_neigh(i); 
    UW((k-1)*3+1:3*k) = UW((k-1)*3+1:3*k)-Kf'*[tau_stokes_nonpx((k-1)*N_f+1+P*N_c:k*N_f+P*N_c); 
        tau_stokes_nonpy((k-1)*N_f+1+P*N_c:k*N_f+P*N_c)];
    Kim = getKmat2D(rimage_k{k},q(k));
    % TODO: don't build K matrix
    UW((k-1)*3+1:3*k) = UW((k-1)*3+1:3*k)-Kim'*[tau_stokes_e_nonpx{k}; tau_stokes_e_nonpy{k}];
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
geom.rcheck = rcheck_b;
u_rhs = matvec_mob_pairprecond_enhanced(tau,geom,basis);

S_0 = getRecompletionFlow(rin,rcheck_b,q,F,T); 
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


%Some visualisation stuff... 
if visualise
% Visualise each component?
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
% sgtitle('Error on boundary mob pair corr','interpreter','latex')

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
    ylabel(c,'Relative residual $\log_{10}$','interpreter','latex','FontSize',14)
    c.FontSize = 14;
    title('Relative residual mob pair corr','interpreter','latex')


    %% Visualise source strengths
    figure()
    semilogy(abs(lambda))
    title('Source strengths mobility with pair corr')
    axis tight
    
end

end


%Old stuff from resistance


function doPairBasisTest(Upf,Ypf,i,p2,q,U,Y,rbase_in_c,rbase_out_f,rpair_fine)
%Test pair basis with a smooth coarse density
warning('Check consistency with the matvec! NOT in use')
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
        step1 = U{i}'*[tau_p1_x;tau_p1_y]; %here I assume x and y follow each other?
        tau_mapped = Y{i}*step1;
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
delta = 0.01;
q = [0; 2+delta; (2+delta)*1i]; %center coordinates
q = [0; 2+delta; 4+4*delta];
%or, instead, three cireles in triangle
% delta = 0.001; 
x = 1+delta/2;
y = sqrt((2+delta)^2-(1+delta/2)^2);
q = [0; 2+delta; x+1i*y];

F = [1 0; 0 0; 0 1]; %forces on the particles
T = [1; 1; 1]; %torques on the particles
rads = [1; 1; 1]; 


%If only two particles
%q = [0; 2+delta; 6];
% F = F(1:2,:); 
% T = T(1:2); 
% rads = [1;1]; 

visualise = 1; 
images = 1; 
delta_pair = 0.5; 
lr= 0; 
%[UW1,lambda_mob,it1,gmres_tol,err1] = solve_2D_mob(q,F,T,rads,images, lr, visualise);

%compare to a solution with image enhancement

[UW2,lambda,it2,gmres_tol,err2] = solve_mob_precond_enhanced(q,F,T,rads,delta_pair,visualise);

str = sprintf('Relative residual with 1-body precond: %1.2e vs 2-body: %1.2e\n Converging in %u resp % u iterations',err1,err2,it1,it2);
disp(str)

alignfigs;

end
