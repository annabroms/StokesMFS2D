function [UW,sol] = solve_mob_1B_enhanced(q,F,T,opt)
%SOLVE_MOB_1B_ENHANCED Solve a 2D Stokes mobility problem with 1-body
%preconditioning and Stokeslet-only ellipse-segment enhancement.
%
% Syntax:
%   [UW,sol] = solve_mob_1B_enhanced(q,F,T,opt)
%
% Inputs:
%   q   - Complex particle centers (P x 1).
%   F   - Px2 net force vectors.
%   T   - Px1 torque vector.
%   opt - Options struct, matching the enhanced solver style. Supports
%         opt.RAM_check for memorygraph-based RAM estimates.
%
% Outputs:
%   UW  - 3P x 1 rigid-body velocity vector [Ux; Uy; W] per particle.
%   sol - Struct with source data, residual information, and
%         sol.ram_estimate.

if nargin==0, test_solve_mob_1B_enhanced; return; end

if nargin < 4 || ~isstruct(opt)
    error('solve_mob_1B_enhanced requires q, F, T, and an options struct opt.');
end

[ram_check,ram_cleanup] = startRamCheck(opt,mfilename); %#ok<NASGU>

q = q(:);
T = T(:);
P = numel(q);
assert(size(F,1)==P,'F must have one row per particle.');
assert(size(F,2)==2,'F must have two columns [Fx, Fy].');
assert(numel(T)==P,'T must have one entry per particle.');

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

opt.visualise_sol = visualise_sol;
opt.gmres_tol = gmres_tol;
opt.gmres_verbose = gmres_verbose;
opt.debug = debug;
opt.maxit = maxit;
opt.get_bndry_field = get_bndry_field;
opt.surface_error_mode = surface_error_mode;
opt.solve_resistance = false;
opt.self_tol = 1e-9; 

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

%% Right-hand side from completion flow
[urhs,~,~] = getRecompletionFlowStokes1BEnhanced(q,F,T,geom);

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
end

ram_check = markRamCheckPhase(ram_check,'precomp_end');

%% Solve
disp(' == Solving... == ');
[tau,it,resvec,real_res] = helsing_gmres( ...
    @(x) matvecStokes1BEnhanced(x,geom,basis), ...
    urhs,2*geom.total_target_count,maxit,gmres_tol,gmres_verbose,geom.rout);
ram_check = markRamCheckPhase(ram_check,'solve_end');

%% Recover source strengths
[lambda_x_raw,lambda_y_raw,lambda_body,lambda_px,lambda_py] = ...
    mapStokes1BEnhancedSources(tau,geom,basis);

UW = zeros(3*P,1);
for k = 1:P
    UW(3*(k-1)+1:3*k) = -bodyMomentsStokes2D(geom.rin_body{k},q(k),lambda_body{k});
end

%% Residual on a separate boundary
if get_bndry_field
    n_bound = 803;
    body_rad = geom.rad;
    rcheck_b = zeros(P*n_bound,1);
    for k = 1:P
        t = linspace(0,2*pi,n_bound+1)';
        t = t(1:end-1);
        inds = (k-1)*n_bound+1:k*n_bound;
        rcheck_b(inds) = q(k) + body_rad*(cos(t)+1i*sin(t));
    end

    geom_check = buildStokes1BCheckGeom(geom,rcheck_b);
    S0 = getRecompletionFlowStokes1BEnhanced(q,F,T,geom_check);
    u_src = getVelocityField(geom.rvec_in,rcheck_b,lambda_px,lambda_py,[],[],[],[],[],[],[],geom.use_fmm);

    u_lhs = zeros(2*numel(rcheck_b),1);
    for k = 1:P
        pts = geom_check.rout_body{k};
        u_k = evaluateRigidMotion2D(q(k),UW(3*(k-1)+1:3*(k-1)+2),UW(3*k),pts);
        u_lhs(geom_check.target_x_ind{k}) = u_k(1:n_bound);
        u_lhs(geom_check.target_y_ind{k}) = u_k(n_bound+1:end);
    end

    u_rhs = u_src - S0;
    diff_vec = u_rhs - u_lhs;
    diff_complex = diff_vec(1:end/2) + 1i*diff_vec(end/2+1:end);
    rhs_complex = S0(1:end/2) + 1i*S0(end/2+1:end);
    boundary_err = abs(diff_complex);
    abs_res = max(boundary_err);
    if any(rhs_complex)
        rel_vec = boundary_err/max(abs(rhs_complex));
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
    title('GMRES convergence, mobility 1B enhanced','interpreter','latex')

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
% lambda_proj is used in the matvec; lambda_raw and lambda_body keep the
% non-projected densities for rigid-motion postprocessing.
sol.lambda = [lambda_px; lambda_py];
sol.lambda_raw = [lambda_x_raw; lambda_y_raw];
sol.lambda_proj = [lambda_px; lambda_py];
sol.lambda_x = lambda_px;
sol.lambda_y = lambda_py;
sol.lambda_x_raw = lambda_x_raw;
sol.lambda_y_raw = lambda_y_raw;
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

function test_solve_mob_1B_enhanced
close all;

delta = 0.001;
q = [0; 2+delta];
F = [1 0; -1 0];
T = [1; -1];
rad = ones(numel(q),1);

opt = get2Dparams(numel(q));
opt.rad = rad(1);
opt.delta_pair = 0.2;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.get_bndry_field = 1;
opt.gmres_tol = 1e-10;
opt.gmres_verbose = 0;
opt.debug = 0;
opt.maxit = 800;
opt.use_fmm = true;

[UW_new,sol_new] = solve_mob_1B_enhanced(q,F,T,opt);
[U_new,W_new] = unpackUW(UW_new);
[FT_back,sol_res] = solve_res_1B_enhanced(q,U_new,W_new,opt);
[F_ref,T_ref] = deal(F,T);
two_way_mob_res = relerr_inf(FT_back,packFT(F_ref,T_ref));
[UW_img,~,it_img,~,rel_img,abs_img] = ...
    solve_mob_1B(q,F,T,rad,1,0,0,opt.gmres_tol,0,'rel',0);

fprintf('Mobility 1B enhanced: rel_res = %.3e, it = %d\n',sol_new.rel_res,sol_new.it);
fprintf('Resistance 1B enhanced (from mobility): rel_res = %.3e, it = %d\n', ...
    sol_res.rel_res,sol_res.it);
fprintf('Legacy image 1B     : rel_res = %.3e, abs_res = %.3e, it = %d\n',rel_img,abs_img,it_img);
fprintf('  physical output relerr (UW) = %.3e\n',relerr(UW_new,UW_img));
fprintf('  two-way error (mob -> res)  = %.3e\n',two_way_mob_res);

end

function e = relerr(a,b)
den = norm(b);
if den == 0
    e = norm(a-b);
else
    e = norm(a-b)/den;
end
end

function e = relerr_inf(a,b)
e = norm(a-b,inf)/max(1,norm(b,inf));
end

function [U,W] = unpackUW(UW)
P = numel(UW)/3;
U = [UW(1:3:end), UW(2:3:end)];
W = UW(3:3:end);
assert(size(U,1)==P,'Unexpected size in unpackUW.');
end

function FT = packFT(F,T)
P = size(F,1);
FT = zeros(3*P,1);
FT(1:3:end) = F(:,1);
FT(2:3:end) = F(:,2);
FT(3:3:end) = T;
end
