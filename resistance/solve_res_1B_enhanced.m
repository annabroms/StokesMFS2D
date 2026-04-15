function [FT,sol] = solve_res_1B_enhanced(q,U,W,opt)
%SOLVE_RES_1B_ENHANCED Solve a 2D Stokes resistance problem with 1-body
%preconditioning and Stokeslet-only ellipse-segment enhancement.
%
% Syntax:
%   [FT,sol] = solve_res_1B_enhanced(q,U,W,opt)
%
% Inputs:
%   q   - Complex particle centers (P x 1).
%   U   - Px2 translational velocities.
%   W   - Px1 angular velocities.
%   opt - Options struct, matching the enhanced solver style. Supports
%         opt.RAM_check for memorygraph-based RAM estimates.
%
% Outputs:
%   FT  - 3P x 1 force/torque vector per particle.
%   sol - Struct with source data, residual information, and
%         sol.ram_estimate.

if nargin==0, test_solve_res_1B_enhanced; return; end

if nargin < 4 || ~isstruct(opt)
    error('solve_res_1B_enhanced requires q, U, W, and an options struct opt.');
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
gmres_verbose = getOptField(opt,'gmres_verbose',0);
debug = logical(getOptField(opt,'debug',false));
maxit = getOptField(opt,'maxit',800);
get_bndry_field = logical(getOptField(opt,'get_bndry_field',true));
get_precomp_time = logical(getOptField(opt,'get_precomp_time',false));
surface_error_mode = lower(char(getOptField(opt,'surface_error_mode','rel')));
if ~any(strcmp(surface_error_mode,{'abs','rel'}))
    error('surface_error_mode must be ''abs'' or ''rel''.');
end
precomp_time = struct('total',nan,'one_body',nan,'two_body_or_peanut',nan);

solver_name = 'res_1B_right';

opt.visualise_sol = visualise_sol;
opt.gmres_tol = gmres_tol;
opt.gmres_verbose = gmres_verbose;
opt.debug = debug;
opt.maxit = maxit;
opt.get_bndry_field = get_bndry_field;
opt.surface_error_mode = surface_error_mode;
opt.solve_resistance = true;

fprintf('==== START: %s ====\n', mfilename);

if get_precomp_time
    one_body_timer = tic;
end
[geom,basis,~,~] = prepareStokes1BEnhanced(q,opt);
if get_precomp_time
    precomp_time.one_body = toc(one_body_timer);
    precomp_time.two_body_or_peanut = 0;
    precomp_time.total = precomp_time.one_body;
end
P = numel(q);

%% Right-hand side from rigid motion on the enhanced collocation grid
fout = zeros(2*geom.total_target_count,1);
for k = 1:P
    u_k = evaluateRigidMotion2D(q(k),U(k,:),W(k),geom.rout_body{k});
    fout(geom.target_x_ind{k}) = u_k(1:geom.target_count(k));
    fout(geom.target_y_ind{k}) = u_k(geom.target_count(k)+1:end);
end

%% Optional matrix inspection
if debug
    ncols = 2*geom.total_target_count;
    CC = zeros(ncols);
    x = zeros(ncols,1);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvecStokes1BEnhanced(x,geom,basis);
    end
    figure();
    imagesc(log10(abs(CC)))
    colorbar
    title([mfilename ': log_{10} |matvec system matrix|'],'interpreter','none')
    cc = skeel(CC);
    fprintf('Estimated condition number of system matrix: %1.3e \n',cc);

    figure();
    [V,D] = eig(CC);
    D = diag(D); 
    plot(real(D),imag(D),'b+')
    title([solver_name ': eigenvalues of matvec system matrix'],'interpreter','none')

    rout = geom.rout;
    num_min = 2;
    [e_sel,ind] = mink(real(D),num_min);
    V_sel = V(:,ind);
    figure()   
    for k = 1:num_min
        subplot(num_min,1,k);
        quiver(real(rout),imag(rout),V_sel(1:end/2,k),V_sel(end/2+1:end,k));
        axis equal
    end
end

ram_check = markRamCheckPhase(ram_check,'precomp_end');

%% Solve
disp(' == Solving... == ');
[tau,it,resvec,real_res] = helsing_gmres( ...
    @(x) matvecStokes1BEnhanced(x,geom,basis), ...
    fout,2*geom.total_target_count,maxit,gmres_tol,gmres_verbose,geom.rout);
ram_check = markRamCheckPhase(ram_check,'solve_end');

%% Recover source strengths
[lambda_x_raw,lambda_y_raw,lambda_body,lambda_px,lambda_py] = ...
    mapStokes1BEnhancedSources(tau,geom,basis);

%% Compute physical output: forces and torques on each body.
FT = zeros(3*P,1);
for k = 1:P
    FT(3*(k-1)+1:3*k) = bodyMomentsStokes2D(geom.rin_body{k},q(k),lambda_body{k});
end

%% Residual on a separate boundary
if get_bndry_field
    n_bound = 803;
    u_true = zeros(2*P*n_bound,1);
    body_rad = geom.rad;
    rcheck_b = zeros(P*n_bound,1);
    for k = 1:P
        t = linspace(0,2*pi,n_bound+1)';
        t = t(1:end-1);
        pts = q(k) + body_rad*(cos(t)+1i*sin(t));
        inds = (k-1)*n_bound+1:k*n_bound;
        rcheck_b(inds) = pts;
        u_k = evaluateRigidMotion2D(q(k),U(k,:),W(k),pts);
        u_true(inds) = u_k(1:n_bound);
        u_true(P*n_bound+inds) = u_k(n_bound+1:end);
    end

    u_src = getVelocityField(geom.rvec_in,rcheck_b,lambda_x_raw,lambda_y_raw,[],[],[],[],[],[],[],geom.use_fmm);
    diff_vec = u_src - u_true;
    diff_complex = diff_vec(1:end/2) + 1i*diff_vec(end/2+1:end);
    true_complex = u_true(1:end/2) + 1i*u_true(end/2+1:end);
    boundary_err = abs(diff_complex);
    abs_res = max(boundary_err);
    if any(true_complex)
        rel_vec = boundary_err/max(abs(true_complex));
    else
        rel_vec = boundary_err;
    end
    rel_res = max(rel_vec);

    fprintf('Relative boundary error: %.3e\n', rel_res);
    fprintf('Absolute boundary error: %.3e\n', abs_res);
else
    rel_res = nan;
    abs_res = nan;
    rel_vec = [];
    boundary_err = [];
    fprintf('Boundary field evaluation skipped (opt.get_bndry_field=0)\n');
end

if visualise_sol
    figure();
    semilogy(resvec);
    axis tight
    grid on
    xlabel('iteration number','interpreter','latex');
    ylabel('Estimated relative residual','interpreter','latex');
    title('GMRES convergence, resistance 1B enhanced','interpreter','latex')

    figure()
    plot([lambda_x_raw; lambda_y_raw]);
    axis tight
    title('Source strength magnitude')


    if get_bndry_field
        if strcmp(surface_error_mode,'rel')
            err_plot = rel_vec;
            err_label = 'relative';
        else
            err_plot = boundary_err;
            err_label = 'absolute';
        end
        figure();
        semilogy(err_plot);
        axis tight
        grid on
        xlabel('check-node index','interpreter','latex');
        ylabel(sprintf('%s boundary error',err_label),'interpreter','latex');
        title(sprintf('%s: %s boundary error',mfilename,err_label),'interpreter','none');
    end
end

sol = struct();
sol.lambda = [lambda_x_raw; lambda_y_raw];
sol.lambda_proj = [lambda_px; lambda_py];
sol.lambda_px = lambda_px;
sol.lambda_py = lambda_py;
sol.lambda_body = lambda_body;
sol.it = it;
sol.gmres_unknowns = 2*geom.total_target_count;
sol.gmres_tol = gmres_tol;
sol.rel_res = rel_res;
sol.abs_res = abs_res;
sol.resvec = resvec;
sol.real_res = real_res;
sol.precomp_time = precomp_time;
sol.ram_estimate = finishRamCheck(ram_check);

end

function test_solve_res_1B_enhanced
close all;

delta = 10; 
q = [0; 2+delta; (2+delta)*1i];
U = [1 0; -1 0; 1 1];
W = [1; -1; 1];

P = 5;
q = grow_cluster(P,delta,2);
U = rand(P,2); 
W = rand(P,1); 
rad = ones(numel(q),1);

N_c = 60; 
N_f = 60; 
opt = get2Dparams(numel(q),N_c,N_f);
opt.rad = rad(1);
opt.delta_pair = 0.2;
opt.visualise_sol = 1;
opt.visualise_grid = 1;
opt.get_bndry_field = 1;
opt.gmres_tol = 1e-10;
opt.gmres_verbose = 0;
opt.debug = 1;
opt.maxit = 800;
opt.use_fmm = true;
%opt.beta = 0.5;

[FT_new,sol_new] = solve_res_1B_enhanced(q,U,W,opt);
debug = 0;
[FT_img,~,it_img,~,rel_img] = ...
    solve_res_1B(q,U,W,rad,1,0,0,opt.gmres_tol,debug,0);

fprintf('Resistance 1B enhanced: rel_res = %.3e, it = %d\n',sol_new.rel_res,sol_new.it);
fprintf('Legacy image 1B      : rel_res = %.3e, it = %d\n',rel_img,it_img);
fprintf('  physical output relerr (FT) = %.3e\n',relerr(FT_new,FT_img));

end

function e = relerr(a,b)
den = norm(b);
if den == 0
    e = norm(a-b);
else
    e = norm(a-b)/den;
end
end
