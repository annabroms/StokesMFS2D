function [UW,sol] = solve_mob_1B_enhanced_left(q,F,T,opt)
%SOLVE_MOB_1B_ENHANCED_LEFT Solve a 2D Stokes mobility problem with
%left-preconditioned 1-body ellipse-segment enhancement.
%
% Syntax:
%   [UW,sol] = solve_mob_1B_enhanced_left(q,F,T,opt)
%
% Inputs:
%   q   - Complex particle centers (P x 1).
%   F   - Px2 net force vectors.
%   T   - Px1 torque vector.
%   opt - Options struct, matching the enhanced solver style.
%
% Outputs:
%   UW  - 3P x 1 rigid-body velocity vector [Ux; Uy; W] per particle.
%   sol - Struct with source data and residual information.

if nargin==0, test_solve_mob_1B_enhanced_left; return; end

if nargin < 4 || ~isstruct(opt)
    error('solve_mob_1B_enhanced_left requires q, F, T, and an options struct opt.');
end

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
surface_error_mode = lower(char(getOptField(opt,'surface_error_mode','rel')));
if ~any(strcmp(surface_error_mode,{'abs','rel'}))
    error('surface_error_mode must be ''abs'' or ''rel''.');
end

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

[geom,basis,~,~] = prepareStokes1BEnhanced(q,opt);
P = numel(q);

%% Right-hand side from completion flow
[urhs,~,~] = getRecompletionFlowStokes1BEnhanced(q,F,T,geom);
rleft = applyStokes1BEnhancedPseudo(urhs,geom,basis);

%% Optional matrix inspection
if debug
    ncols = 2*geom.total_source_count;
    CC = zeros(ncols);
    x = zeros(ncols,1);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvecStokes1BEnhancedMobilityLeft(x,geom,basis);
    end
    figure();
    imagesc(log10(abs(CC)))
    colorbar
    title([mfilename ': log_{10} |matvec system matrix|'],'interpreter','none')
    cc = skeel(CC);
    fprintf('Estimated condition number of system matrix: %1.3e \n',cc);

    e = eig(CC);
    figure()
    scatter(real(e),imag(e),'+');
    title('Eigenvalues of left-preconditioned mobility system matrix')
end

%% Solve
disp(' == Solving... == ');
[lambda,it,resvec,real_res] = helsing_gmres( ...
    @(x) matvecStokes1BEnhancedMobilityLeft(x,geom,basis), ...
    rleft,2*geom.total_source_count,maxit,gmres_tol,gmres_verbose,geom.rvec_in);

[lambda_x_raw,lambda_y_raw,lambda_body,lambda_px,lambda_py] = ...
    unpackProjectedMobilityLambda(lambda,geom);

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
    title('GMRES convergence, mobility 1B enhanced left','interpreter','latex')

    figure()
    plot([lambda_x_raw; lambda_y_raw]);
    axis tight
    title('Raw source strength magnitude')

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
% The left solve uses the raw source density as the GMRES unknown while the
% projected density is the one used in the physical Stokeslet field.
sol.lambda = lambda;
sol.lambda_raw = lambda;
sol.lambda_proj = [lambda_px; lambda_py];
sol.lambda_x = lambda_x_raw;
sol.lambda_y = lambda_y_raw;
sol.lambda_x_raw = lambda_x_raw;
sol.lambda_y_raw = lambda_y_raw;
sol.lambda_px = lambda_px;
sol.lambda_py = lambda_py;
sol.lambda_body = lambda_body;
sol.it = it;
sol.gmres_tol = gmres_tol;
sol.rel_res = rel_res;
sol.abs_res = abs_res;
sol.resvec = resvec;
sol.real_res = real_res;

end

function test_solve_mob_1B_enhanced_left
close all;

delta = 0.001;
q = [0; 2+delta];
P = 5; 
warning(['delta small and P e.g 10 will require a lot of iterations, ' ...
    'particularly for resistance!'])
q = grow_cluster(P,delta,2);
F = rand(P,2);
T = rand(P,1);
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
opt.maxit = 1000; %max number of iterations
opt.use_fmm = true;

[field_err,identity_err] = run_left_preconditioner_sanity(opt);
[UW_left,sol_left] = solve_mob_1B_enhanced_left(q,F,T,opt);
[UW_right,sol_right] = solve_mob_1B_enhanced(q,F,T,opt);
left_right_err = relerr(UW_left,UW_right);
FT_ref = packFT(F,T);

[U_left,W_left] = unpackUW(UW_left);
[FT_back_left,sol_res_left] = solve_res_1B_enhanced_left(q,U_left,W_left,opt);
two_way_mob_res_left = relerr_inf(FT_back_left,FT_ref);

[U_right,W_right] = unpackUW(UW_right);
[FT_back_right,sol_res_right] = solve_res_1B_enhanced_left(q,U_right,W_right,opt);
two_way_mob_res_right = relerr_inf(FT_back_right,FT_ref);

fprintf('Left mobility pseudo sanity: rel_field = %.3e, rel_id = %.3e\n', ...
    field_err,identity_err);
fprintf('Mobility 1B left         : rel_res = %.3e, it = %d\n', ...
    sol_left.rel_res,sol_left.it);
fprintf('Mobility 1B right        : rel_res = %.3e, it = %d\n', ...
    sol_right.rel_res,sol_right.it);
fprintf('Resistance 1B left (from left mob)  : rel_res = %.3e, it = %d\n', ...
    sol_res_left.rel_res,sol_res_left.it);
fprintf('Resistance 1B left (from right mob) : rel_res = %.3e, it = %d\n', ...
    sol_res_right.rel_res,sol_res_right.it);
fprintf('  physical output relerr (left vs right) = %.3e\n', ...
    left_right_err);
fprintf('  two-way error   (left mob -> left res) = %.3e\n', ...
    two_way_mob_res_left);
fprintf('  two-way error  (right mob -> left res) = %.3e\n', ...
    two_way_mob_res_right);


end

function [field_err,identity_err] = run_left_preconditioner_sanity(opt)
opt_local = opt;
opt_local.visualise_sol = 0;
opt_local.visualise_grid = 0;
opt_local.get_bndry_field = 0;
opt_local.debug = 0;
opt_local.use_fmm = false;
opt_local.solve_resistance = false;

[geom,basis] = prepareStokes1BEnhanced(0,opt_local);
rng(1);
u_rand = randn(2*geom.total_target_count,1);
lambda_ref = basis.Y{1}*(basis.U{1}*u_rand);
u_self = buildLeftMobilityBoundaryField(lambda_ref,geom);
lambda_rec = applyStokes1BEnhancedPseudo(u_self,geom,basis);

[~,~,lambda_body] = unpackStokes1BEnhancedLambda(lambda_ref,geom);
mom = bodyMomentsStokes2D(geom.rin_body{1},geom.q(1),lambda_body{1});
u_direct = basis.Aii{1}*lambda_body{1} + ...
    evaluateRigidMotion2D(geom.q(1),mom(1:2).',mom(3),geom.rout_body{1});

field_err = relerr(u_self,u_direct);
identity_err = relerr(lambda_rec,lambda_ref);
end

function res = matvecStokes1BEnhancedMobilityLeft(lambda,geom,basis)
u = buildLeftMobilityBoundaryField(lambda,geom);
res = applyStokes1BEnhancedPseudo(u,geom,basis);
end

function u = buildLeftMobilityBoundaryField(lambda,geom)
[~,~,lambda_body,lambda_px,lambda_py] = unpackProjectedMobilityLambda(lambda,geom);

u = getVelocityField(geom.rvec_in,geom.rout,lambda_px,lambda_py, ...
    [],[],[],[],[],[],[],geom.use_fmm);

for k = 1:numel(geom.q)
    mom = bodyMomentsStokes2D(geom.rin_body{k},geom.q(k),lambda_body{k});
    u_k = evaluateRigidMotion2D(geom.q(k),mom(1:2).',mom(3),geom.rout_body{k});
    nt = geom.target_count(k);
    u(geom.target_x_ind{k}) = u(geom.target_x_ind{k}) + u_k(1:nt);
    u(geom.target_y_ind{k}) = u(geom.target_y_ind{k}) + u_k(nt+1:end);
end
end

function [lambda_x_raw,lambda_y_raw,lambda_body,lambda_px,lambda_py] = ...
    unpackProjectedMobilityLambda(lambda,geom)
[lambda_x_raw,lambda_y_raw,lambda_body] = unpackStokes1BEnhancedLambda(lambda,geom);

lambda_px = zeros(geom.total_source_count,1);
lambda_py = zeros(geom.total_source_count,1);

for k = 1:numel(geom.q)
    ns = geom.source_count(k);
    sx = geom.source_x_ind{k};
    lambda_k_proj = projectOutRigid2D(lambda_body{k},geom.rin_body{k},geom.q(k));
    lambda_px(sx) = lambda_k_proj(1:ns);
    lambda_py(sx) = lambda_k_proj(ns+1:end);
end
end

function [lambda_x,lambda_y,lambda_body] = unpackStokes1BEnhancedLambda(lambda,geom)
lambda = lambda(:);
ns_tot = geom.total_source_count;
if numel(lambda) ~= 2*ns_tot
    error('unpackStokes1BEnhancedLambda:badInputSize', ...
        'lambda must have length 2*geom.total_source_count.');
end

lambda_x = lambda(1:ns_tot);
lambda_y = lambda(ns_tot+1:end);
lambda_body = cell(numel(geom.q),1);
for k = 1:numel(geom.q)
    sx = geom.source_x_ind{k};
    lambda_body{k} = [lambda_x(sx); lambda_y(sx)];
end
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
