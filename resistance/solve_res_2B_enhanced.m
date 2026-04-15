function [FT,sol] = solve_res_2B_enhanced(q,U,W,opt)
%SOLVE_RES_PRECOND_ENHANCED Solves a 2D Stokes resistance problem with circular
%particles using MFS with 2-body preconditioning. To resolve
% challenging close interactions, a fine 2-body BVP is solved for fine
% source strengths on each close pair of particles. The fine sources,
% including shielding of image points, correct the representation
% obtained from a coarse grid, effectively preconditioning the system.
%
% Syntax:
%   [FT,sol] = solve_res_2B_enhanced(q,U,W,opt)
%
% Inputs:
%   q          - Vector of length P, complex-valued center coordinates for the particles
%   U          - Px2 matrix of translational velocities (columns: x and y components)
%   W          - Px1 column vector of angular velocities
%   opt        - Options struct. Common fields:
%                rad, delta_pair, lr, visualise_sol, gmres_tol,opt.gmres_verbose,
%                use_fmm, N_c, N_f, a_c, a_f, tol_c, RAM_check.
%       debug    build/plot/investigate system matrix corresponding to
%                matvec.
%
% Outputs:
%   FT         - 3P×1 vector of computed net forces and torques 
%   sol        - Struct with fields:
%                lambda, it, gmres_tol, rel_res, resvec, ram_estimate.
%
% Description:
%   This function applies a 2-body preconditioner (using pair corrections via local fine BVPs)
%   to improve solution accuracy for near-contact particle configurations.
%   The FMM is used for Stokeslet evaluation.
%
% Notes:
%   - This is an MFS generalisation of the Cheng-Greengard idea (1998/2000)
%   - Intended to match solve_mob_2B_images, but applied to a resistance problem rather than mobility.
%
% See also:
%   solve_res_1B              - 1-body preconditioned resistance solver
%   solve_res_peanut_enhanced - 2-body preconditioner with peanut
%                               compression (in this version of the algorithm, only coarse sources
%                               are needed at the solve stage).
%   solve_mob_2B_enhanced      - 2-body preconditioned mobility solver
%   matvec_res_2B_enhanced     - matrix-vector action used by GMRES
%
% To test: Call without arguments.
%
% Anna Broms, Mar 2026

if nargin==0, test_solve_res; 
    return; end

if nargin < 4 || ~isstruct(opt)
    error('solve_res_2B_enhanced requires q, U, W, and an options struct opt.');
end

[ram_check,ram_cleanup] = startRamCheck(opt,mfilename); %#ok<NASGU>

q = q(:);
W = W(:);
P = numel(q);
assert(size(U,1)==P,'U must have one row per particle.');
assert(size(U,2)==2,'U must have two columns [Ux, Uy].');
assert(numel(W)==P,'W must have one entry per particle.');

visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = getOptField(opt,'gmres_tol',1e-10);
debug = logical(getOptField(opt,'debug',false));
gmres_verbose = getOptField(opt,'gmres_verbose',0);
maxit = getOptField(opt,'maxit',800);
use_fmm = logical(getOptField(opt,'use_fmm',true));
lr = getOptField(opt,'lr',0);
get_precomp_time = logical(getOptField(opt,'get_precomp_time',false));
rad = ones(P,1); %assume for now all are identical
precomp_time = struct('total',nan,'one_body',nan,'pair_setup',nan, ...
    'pair_basis',nan,'two_body_or_peanut',nan);

%% SET PARAMS
if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

opt.P = P;
opt.gmres_verbose = gmres_verbose;
opt.use_fmm = use_fmm;


%Play with N_c, N_f, a (a_f). 
N_c = getOptField(opt,'N_c',60);
N_f = getOptField(opt,'N_f',60);
a_c = getOptField(opt,'a_c',1.2);
a_f = getOptField(opt,'a_f',1.2);

% Determine separation to proxy surface.
tol_c = getOptField(opt,'tol_c',1e-12);

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rp_c = getOptField(opt,'Rp_c',max([1-sep_c,0.01]));
Rp_f = getOptField(opt,'Rp_f',max([1-sep_f,0.01]));

 
%accumulation point, given Rp and delta. Closed formula from fixed point of reflection formula
accstop = (1-Rp_c)^2/Rp_c;  

delta_pair = getOptField(opt,'delta_pair',accstop);
 
opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c;
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 

opt.P = P; 
opt.N_peanut = 0; 
opt.precomp = getOptField(opt,'precomp',1);
opt.pc = 1; %prepare grid to do pair corrections
opt.delta_pair = delta_pair; 
opt.lr = lr; %using long range preconditioning? 
opt.Nclust = getOptField(opt,'Nclust',200);
opt.beta = getOptField(opt,'beta',0.3);


%% CREATE GRID
%Outer basic grid
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';
rbase_out_c = rad(1)*cos(tout_c)+1i*rad(1)*sin(tout_c);
%Inner basic grid
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rp_c*cos(tin)+Rp_c*1i*sin(tin); 
tin_f = linspace(0,2*pi,N_f+1);
tin_f = tin_f(1:end-1)'; 
rbase_in_f = Rp_f*cos(tin_f)+Rp_f*1i*sin(tin_f); 

%Construct image grid
%will return only the basic outer grid, else refined outer grid 

rout = []; 
for k = 1:P
    rout = [rout; rbase_out_c+q(k)];
end
 
if get_precomp_time
    pair_setup_timer = tic;
end
[~, ~, ~, rimage_vec, refine,pairs] = getEnhancedGrid(q, opt);
if get_precomp_time
    precomp_time.pair_setup = toc(pair_setup_timer);
end

%Get pair basis
if get_precomp_time
    pair_basis_timer = tic;
end
[Upf,Ypf,~,~,~,Cmap_FU,pair_cache] = ...
    getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,[],rbase_out_c);
if get_precomp_time
    precomp_time.pair_basis = toc(pair_basis_timer);
    precomp_time.two_body_or_peanut = precomp_time.pair_setup + precomp_time.pair_basis;
end


%Get one-body pseduo inverse blocks -- enough to do this for single body (everybody has the same coarse grid).
if get_precomp_time
    one_body_timer = tic;
end
[UU,YY] = getSelfPseudo(1,rbase_in_c,rbase_out_c);
if get_precomp_time
    precomp_time.one_body = toc(one_body_timer);
    precomp_time.total = precomp_time.one_body + precomp_time.two_body_or_peanut;
end

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.rimage_vec = rimage_vec;
geom.opt = opt;
geom.rvec_out = rout;
geom.q = q;
geom.pairs = pairs;

basis = struct();
basis.U = UU;
basis.Y = YY;
basis.Upf = Upf;
basis.Ypf = Ypf;
basis.Cmap_FU = Cmap_FU;
basis.pair_cache = pair_cache;


%Visualise 1-body and pair-basis

%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,[],[],N_c, N_f,a_c,a_f,rad)

%% Construct rhs

%Set boundary condition for MFS. Evaluate known velocities on the coarse
%grid of collocation points. 
foutx = [];
fouty = [];
for k = 1:P
   rhs_x = @(x) U(k,1)-W(k)*(imag(x)-imag(q(k)));
   rhs_y = @(x) U(k,2)+W(k)*(real(x-q(k)));
   foutx = [foutx; rhs_x(rout((k-1)*ceil(a_c*N_c)+1:k*ceil(a_c*N_c)))];
   fouty = [fouty; rhs_y(rout((k-1)*ceil(a_c*N_c)+1:k*ceil(a_c*N_c)))];  
end

fout = [foutx; fouty];


%% Construct check boundaries
% Create new grid points, for which the accuracy of the solution is
% to be evaluated. 
rcheck_b = [];
n_bound = 2000; %check in a large number of points
t = linspace(0,2*pi,n_bound)';
for k = 1:P
    rcheck_b = [rcheck_b; q(k)+rad(k)*(cos(t)+1i*sin(t))];
end

%% SOLVE SYSTEM
% Build the matrix to inspect conditioning/eigenvalues if requested.
if debug
    x = zeros(2*length(rout),1);
    ncols = 2*length(rout);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        uu = matvec_res_2B_enhanced(x,geom,basis);
        CC(:,k) = uu;
    end
    figure(13);
    clf; 
    imagesc(log10(abs(CC)))
    colorbar
    title([solver_name ': log_{10} |matvec system matrix|'],'interpreter','none')
    cc = skeel(CC);
    fprintf('Estimated condition number of system matrix: %1.3e \n',cc);

    figure(5)
    [V,D] = eig(CC);
    D = diag(D); 
    plot(real(D),imag(D),'ro')
    title([solver_name ': eigenvalues of matvec system matrix'],'interpreter','none')

    [s,I] = mink(abs(D),3);
    Vsmall = V(:,I).*s';

    
    Mc = round(a_c*N_c);
    t = linspace(0,2*pi,Mc+1);
    t = t(1:end-1)'; 
    
    %% Visualise eigvecs of system matrix
    eigvec_vis = 0; 
    if eigvec_vis
        for i = 1:3
            for k = 1:P  
                Vpx = Vsmall((k-1)*Mc+1:k*Mc,i);
                Vpy = Vsmall((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
                Vp = abs(Vpx+1i*Vpy);
                figure(31)
                subplot(1,3,i)
                scatter3(real(q(k))+sin(t),imag(q(k))+cos(t),Vp,40,Vp,'filled');
                hold on
                colorbar
                view(0,90)
    
                figure(32)
                subplot(1,3,i)
                plot(t,Vp)
                hold on
    
    
                figure(27)
                subplot(1,3,i)
                plot(t,Vpx)
                hold on
    
                figure(28)
                subplot(1,3,i)
                plot(t,Vpy)
                hold on
    
                figure(29)
                subplot(1,3,i)
                scatter3(real(q(k))+sin(t),imag(q(k))+cos(t),Vpx,40,Vpx,'filled');
                hold on
                colorbar
                view(0,90)
    
                figure(30)
                subplot(1,3,i)
                scatter3(real(q(k))+sin(t),imag(q(k))+cos(t),Vpy,40,Vpy,'filled');
                hold on
                colorbar
                view(0,90)
    
    
    
            end

        end
    end
end
%% Experiment with left preconditioner based on deflation

if lr
    warning('Long range preconditioning not tested')
    rin_c = []; 
    for k = 1:P
        rin_c = [rin_c; rbase_in_c+q(k)];
    end
    [Sinv,Nx,Ny,Mx,Z,Y,db] = get_long_range_precond(q,rin_c,rout,opt);
    %tau_coarse1 = AN*Rinv*(AM'*fout); same thing
    opt.db = db;
    tau_coarse = getCoarseSource(fout,Sinv,Nx,Ny,Mx,Z,Y,db,P,opt.N_c,opt.a_c);
end

ram_check = markRamCheckPhase(ram_check,'precomp_end');

disp(' == Solving... == ');
[tau,it,resvec,real_res] = helsing_gmres(@(x) matvec_res_2B_enhanced(x,geom,basis),fout,2*size(rout,1),maxit,gmres_tol,opt.gmres_verbose,rout);
ram_check = markRamCheckPhase(ram_check,'solve_end');

plot_gmres = visualise_sol;

%With Krylov precond, do something like
%[tau, e2, precond] = precond_gmres(@(x) matvec_res_2B_images(x,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rvec_out,q,UU,YY,pairs,Upf,Ypf,s), fout, zeros(2*size(rvec_out,1),1), 2*size(rvec_out,1), gmres_tol, precond,debug);
%it = length(e2); 

if plot_gmres
      figure()
%     semilogy(e2);
      semilogy(resvec);
      xlabel('iteration number','interpreter','latex');
      ylabel('Estimated relative residual');
      axis tight
      grid on
      title('Convergence resistance with pair corr','interpreter','latex')

      %what's the resiudal?
    %  u = matvec_res_2B_images(tau,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,Upf,Ypf,s);
      u = matvec_res_2B_enhanced(tau,geom,basis);
%       figure()
%       semilogy(abs(u-fout))
%       title('Residual in solution')
end 


disp(' == Postprocessing == ');
%% POSTPROCESS
[rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y] = getPairTransformationStokes(tau,geom,basis);

lambda =  [tau_stokes_x; tau_stokes_y];

%% Compute net forces and torque
% FT = [F_x; F_y; T_z] per particle, assembled from coarse, fine and
% pair-specific fine source strengths.
FT = zeros(3*P,1);

% Coarse source contribution.
Kc = getKmat2D(rbase_in_c,0);
for k = 1:P
    seg = (k-1)*3+1:3*k;
    idx = coarse_ind{k};
    FT(seg) = FT(seg) + Kc'*[tau_stokes_x(idx); tau_stokes_y(idx)];
end

if logical(getOptField(opt,'cmap',false)) && (~isempty(basis.Cmap_FU) || basis.pair_cache.enabled)
    for row = 1:size(pairs,1)
        i1 = pairs(row,1);
        i2 = pairs(row,2);
        rhs_pair = [tau_stokes_x(coarse_ind{i1}); tau_stokes_x(coarse_ind{i2}); ...
                    tau_stokes_y(coarse_ind{i1}); tau_stokes_y(coarse_ind{i2})];
        if isfield(basis,'pair_cache') && basis.pair_cache.enabled
            pair = getStokesPairInstance(basis.pair_cache,row);
            rhs_pair = rotatePairOrderedStokesData(rhs_pair,opt.N_c,pair.meta.phase_c,conj(pair.meta.rot));
            ft_pair = pair.group.Cmap_FU*rhs_pair;
            fi = pair.meta.rot*(ft_pair(1) + 1i*ft_pair(2));
            fj = pair.meta.rot*(ft_pair(4) + 1i*ft_pair(5));
            ft_pair = [real(fi); imag(fi); ft_pair(3); ...
                       real(fj); imag(fj); ft_pair(6)];
        else
            ft_pair = basis.Cmap_FU{i1,i2}*rhs_pair;
        end
        FT((i1-1)*3+1:3*i1) = FT((i1-1)*3+1:3*i1) + ft_pair(1:3);
        FT((i2-1)*3+1:3*i2) = FT((i2-1)*3+1:3*i2) + ft_pair(4:6);
    end
else
    % Fine and image-source contributions are stored particle-by-particle in
    % the same order as in getPairTransformationStokes.
    has_neigh = sort(unique(pairs(:)));
    Kf = getKmat2D(rbase_in_f,0);
    offset = P*N_c;

    % Number of (extra) sources per particle for index bookkeeping.
    n_im = zeros(P,1);
    for row = 1:size(pairs,1)
        i1 = pairs(row,1);
        i2 = pairs(row,2);
        n_im(i1) = n_im(i1) + length(rimage_vec{i1,i2});
        n_im(i2) = n_im(i2) + length(rimage_vec{i2,i1});
    end

    for i = 1:length(has_neigh)
        k = has_neigh(i);
        seg = (k-1)*3+1:3*k;

        fine_idx = offset + (1:N_f);
        offset = offset + N_f;
        FT(seg) = FT(seg) + Kf'*[tau_stokes_x(fine_idx); tau_stokes_y(fine_idx)];

        nimk = n_im(k);
        if nimk > 0
            extra_idx = offset + (1:nimk);
            offset = offset + nimk;
            Kim = getKmat2D(rvec_in(extra_idx),q(k));
            FT(seg) = FT(seg) + Kim'*[tau_stokes_x(extra_idx); tau_stokes_y(extra_idx)];
        end
    end

    if offset ~= length(tau_stokes_x)
        error('Force postprocessing indexing mismatch in solve_res_2B_enhanced.');
    end
end


%extract force and torque separately
for k = 1:P
    F(k) = FT((k-1)*3+1) + 1i*FT((k-1)*3+2);
    T(k) = FT(3*k); 
end


%% Do the evaluation of the flow in check points , FMM is applied.

ftest_b = getVelocityField(rvec_in,rcheck_b,...
    tau_stokes_x,tau_stokes_y);

%% Compute residual at boundary
fbound_x = ftest_b(1:length(rcheck_b));
fbound_y = ftest_b(length(rcheck_b)+1:end);

fb_x = [];
fb_y = [];
 
 
for k = 1:P
    rhs_f = @(x) [U(k,1)-W(k)*(imag(x)-imag(q(k))); U(k,2)+W(k)*(real(x-q(k)))]; 
    fb_true = rhs_f(rcheck_b(n_bound*(k-1)+1:n_bound*k));
    fb_x = [fb_x; fb_true(1:n_bound)];
    fb_y = [fb_y; fb_true(n_bound+1:end)];   
end

rel_res = max(sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./max(sqrt(fb_x.^2+fb_y.^2));
fprintf('Max surface rel residual at new nodes is %.3e\n', rel_res);

%Some visualisation stuff... 
if visualise_sol

    figure()
    plot(ftest_b)
    hold on
    plot([fb_x; fb_y])
    title('Velocity field rhs and lhs, res pair corr')

    %Visualise component by component
%     figure(9)
%     subplot(2,2,1)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fb_x-fbound_x)),30,log10(abs(fb_x-fbound_x)),'filled')
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
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fb_y-fbound_y)),30,log10(abs(fb_y-fbound_y)),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('error in y velocity')
% 
%     % Visualise the actual velocity
%     subplot(2,2,3)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fb_x)),30,log10(abs(fb_x)),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('x velocity MFS')
% 
%     subplot(2,2,4)
%     scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(fbound_x)),30,log10(abs(fbound_x)),'filled')
%     colorbar
%     axis equal
%     view(0,90)
%     grid off
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     title('x velocity reference')
% 
% 
%     sgtitle('Error on boundary','interpreter','latex')

    % Visualise the error for x and y together
    figure(10)
    m = max(sqrt(fb_x.^2+fb_y.^2));
    scatter3(real(rcheck_b),imag(rcheck_b),log10((sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./m),...
        30,log10((sqrt((fb_x-fbound_x).^2+(fb_y-fbound_y).^2))./m),'filled')
    c = colorbar;
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    ylabel(c,'Relative error in $u$ on boundary (log10)','Interpreter','latex')
    set(c,'TickLabelInterpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    set(gcf,'color','w');
    axis off
    title('Resistance pair corr')

    %% Visualise density
    figure()
    semilogy(abs(lambda))
    ylabel('Magnitude of MFS coefficients','interpreter','latex')
    title('Pair corr resistance')
    axis tight

    % lambda_x = lambda(1:length(rin_f));
    % lambda_y = lambda(length(rin_f)+1:2*length(rin_f));
    % 
    % figure(12)
    % p = length(rin_f)/length(q);
    % 
    % clf; 
    % lambda_tot = vecnorm([lambda_y lambda_x],2,2);
    % for k = 1:length(q)
    %     scatter(real(rin_f((k-1)*p+1:k*p)),imag(rin_f((k-1)*p+1:k*p)),...
    %         30,log10(abs(lambda_tot((k-1)*p+1:k*p))),'filled');
    %     hold on
    % end
    % %Visualise extra singularities at image points
    % 
    % c = colorbar;
    % axis equal
    % view(0,90)
    % grid off
    % set(gca,'xtick',[])
    % set(gca,'ytick',[])
    % hold on
    % plot(real(rout),imag(rout),'k.')
    % ylabel(c,'Magnitude of MFS coefficients. log10','interpreter','latex')
    % set(c,'TickLabelInterpreter','latex')
    % set(gca,'TickLabelInterpreter','latex')
    % axis off
    % set(gcf,'color','w');
end

sol = struct();
sol.lambda = lambda;
sol.it = it;
sol.gmres_tol = gmres_tol;
sol.rel_res = rel_res;
sol.resvec = resvec;
sol.precomp_time = precomp_time;
sol.pair_precomp_stats = pair_cache.stats;
sol.ram_estimate = finishRamCheck(ram_check);

end


function doPairBasisTest(Upf,Ypf,i,p2,q,U,Y,rbase_in_c,rbase_out_f,rpair_fine)
%Test pair basis with a smooth coarse density
warning('Check consistency with the matvec!')
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
        Nother = stokSLPmat(rbase_in_c+q(i),rbase_out_f+q(p2),mu);
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
        Nother = stokSLPmat(rbase_in_c+q(p2),rbase_out_f+q(i),mu);
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

        Npair = stokSLPmat(rpair_fine,rcheck,mu);
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
        Npair = stokSLPmat(rpair_fine,[rbase_out_f+q(i); rbase_out_f+q(p2)],mu);
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

function test_solve_res

close all; 
delta = 0.01; 
q = [0; 2+delta; (2+delta)*1i]; %center coordinates

%or, instead, three cireles in triangle - this DOES NOT converge to a good
%sol! 
% x = 1+delta/2;
% y = sqrt((2+delta)^2-(1+delta/2)^2);
% q = [0; 2+delta; x+1i*y];

U = [1 0; 0 0; 0 1]; %translational velocities
W = [1; 1; 1]; %angular velocities

%or, a pair of circles
% q = [0; 2+delta]; 
% U = [1 0; 0 0];
% W = [0; 0]; 

rad = [1; 1; 1]; 
visualise = 1; 
images = 1; %only relevant for 1-body precond
delta_pair = 0.2; 
test = 1; % 1 or 2
P = length(q); 
%% compare to a solution with 1 body precond only
if test == 1
    %[FT,lambda,it1,gmres_tol,err1] = solve_res_1B(q,U,W,rad,images, visualise);
    rng(9);
    P = 5;
    delta = 0.05; %P = 5
    delta = 1e-2;
    %delta = 1;
    q = [0; 2+delta]*1i;
    %q = grow_cluster(P,delta,2);
    % side = 2 + delta;               % neighbor center distance
    % R = side / (2*sin(pi/P));         % ring radius
    % q = R * exp(1i * (0:P-1).' * (2*pi/P));
    % q = [0; 2+delta; 4+2*delta]*1i;
    P = length(q);
   % q = [q; -6+1.5i; -2-4i]; P = P+2;
    %q = q([1,2,4],:); P = 3; 
    U = rand(P,2); W = rand(P,1); rad = ones(P,1);
    lr = 20; 
    lr = 0; 
    images = 1; 

   % [FT1,lambda,it1,gmres_tol,err1] = solve_res_1B(q,U,W,rad,images, lr,visualise);
    gmres_tol = 1e-7;
    debug = 1; 
    opt = get2Dparams(P);
    opt.rad = rad;
    opt.delta_pair = delta_pair;
    opt.lr = lr;
    opt.visualise_sol = visualise;
    opt.gmres_tol = gmres_tol;
    opt.debug = debug;
    [FT2,sol2] = solve_res_2B_enhanced(q,U,W,opt);
    [FT3,lambda,it3,gmres_tol,err3] = solve_res_2B_images(q,U,W,rad,delta_pair,lr,visualise);

    rel_FT2_vs_FT1 = norm(FT2-FT1,inf)/max(1,norm(FT1,inf));
    rel_FT3_vs_FT1 = norm(FT3-FT1,inf)/max(1,norm(FT1,inf));
    rel_FT3_vs_FT2 = norm(FT3-FT2,inf)/max(1,norm(FT2,inf));

    fprintf('Relative residuals: 1B=%1.2e, 2B-enhanced=%1.2e, 2B-images=%1.2e\n',err1,sol2.rel_res,err3);
    fprintf('GMRES iterations  : 1B=%u, 2B-enhanced=%u, 2B-images=%u\n',it1,sol2.it,it3);
    fprintf('Relative force / torque errors: ||FT2-FT1||/||FT1|| = %1.2e\n',rel_FT2_vs_FT1);
    fprintf('                    ||FT3-FT1||/||FT1|| = %1.2e\n',rel_FT3_vs_FT1);
    fprintf('                    ||FT3-FT2||/||FT2|| = %1.2e\n',rel_FT3_vs_FT2);


    
    alignfigs;
else
    
    %% determine 2-way error (solve resistance followed by mob)
    lr = 0;
    [FT,lambda,it1,gmres_res, err_res] = solve_res_2B_images(q,U,W,rad,delta_pair,lr,visualise);
    F = [FT(1:3:end) FT(2:3:end)];
    T = FT(3:3:end); 
    [UW,lambdahat,it1,gmres_mob, err_mob] = solve_mob_2B_images(q,F,T,rad,delta_pair,visualise);
    Ures = [U W]';
    str = sprintf('Two way error is %1.3e, with resistance residual %1.3e and mobility residual %1.3e',norm(Ures(:)-UW),err_res,err_mob)
    alignfigs;
end
end
