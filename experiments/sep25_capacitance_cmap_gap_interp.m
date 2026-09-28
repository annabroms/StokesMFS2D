%% SEP25_CAPACITANCE_CMAP_GAP_INTERP Interpolation of the canonical aligned
% Laplace capacitance pair-correction matrix C(delta) in the gap.
%
% Scope: this experiment is restricted to the scalar Laplace capacitance
% problem only (Stokes is explicitly out of scope here). "A(delta)" in the
% task description is the coarse-to-coarse pair map Cmap{1,2} built by
% getPairBasisLaplace with opt.cmap=1: the fixed-size (2*N_c x 2*N_c) matrix
% actually inserted into the global coarse system after the fine two-body
% solve and peanut compression (not either intermediate pseudoinverse
% factor). It is referred to as C(delta) below, following the naming in
% sep24_capacitance_cmap_delta_smoothness.m and
% sep25_capacitance_cmap_rank_truncation.m, which this script builds on and
% reuses directly (buildCanonicalPairCmap.m wraps exactly the same calls).
%
% Canonical aligned pair: two identical discs of radius R, centers on the
% real axis at 0 and 2R+delta, aligned surface/proxy grids, fixed body and
% node ordering. Rotations and non-aligned grids are out of scope.
%
% Two dimensionless gap parametrisations are compared throughout:
%   s     = log(delta/R)
%   alpha = acosh(1 + delta/(2*R))
%
% The fine-pair and peanut pseudoinverses use Tikhonov regularisation, so
% singular values crossing a hard TSVD threshold no longer create rank
% switches. The gap range is split for comparison at
%
%   delta_star = (R-Rp_f)^2/Rp_f,
%
% which was the ellipse-segment switch for the old gap-following
% discretisation.  Here ellipse_constant=true freezes the 150-node
% candidate discretisation at delta_min, so the ellipse sources remain
% active and their retained set is constant throughout both comparison
% ranges.  Exactly two base ranges are still used for each parametrisation.
%
% This script does not modify any production solver file. It only calls
% existing routines (getEnhancedGrid, getPairBasisLaplace, lapSLPmat,
% createPeanut) through the small read-only helper buildCanonicalPairCmap.m.
%
% Anna Broms, Sep 25, 2026

close all;
clearvars;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;
addpath(fullfile(repo_root,'experiments','cap_interp'));

fprintf('=== sep25_capacitance_cmap_gap_interp ===\n\n');

%% ------------------------------------------------------------------
%  Configuration (all tunable, all visible here)
%  ------------------------------------------------------------------
R = 2;                      % disc radius
P = 2;                      % two bodies (canonical aligned pair)

opt = getLaplace2Dparams(P,R);
opt.cmap = 1;                           % build the coarse-to-coarse pair map
opt.compress_cmap = false;              % this separate operation remains TSVD
opt.reuse_pair_basis_by_sep = false;    % build Cmap directly for each gap
opt.parallel_precomp = false;
opt.show_counter = 0;
opt.delta_pair = 0.2*R;                 % production pair-correction cutoff
opt.use_tikhonov = true;
opt.tikhonov_tol = 1e-11;               % lambda/sigma_max; Sep 25 sweep choice
opt.ellipse_constant = true;             % freeze ellipse nodes at delta_min
opt.Nclust = 150;                        % candidate nodes on each segment

% Near-pair interval used by the production method (delta_min,delta_max]:
delta_min = 1e-3*R;                % smallest gap ever requested in production (opt.smallest_delta)
delta_max = 0.999*opt.delta_pair;  % largest gap for which pair corrections are applied
                                    % (getEnhancedGrid requires gap < opt.delta_pair strictly)
opt.smallest_delta = delta_min;

% Reference used in C(delta)-C_ref for the spatial reduction.  Supported
% values are 'minimum', 'midpoint' (arithmetic in physical delta), and
% 'maximum'.  This choice does not change direct barycentric interpolation.
reference_location = 'midpoint';

% Nested Chebyshev-Lobatto refinement levels (must be 2^k+1 to nest exactly).
q_levels = [9 17 33 65];

% Tunable truncation levels (all configurable from here, per the task):
rank_trunc_tol   = 1e-8;   % spatial-rank truncation level for section 4/5 (sep25-style)
matrix_action_tol = 1e-7;  % default target matrix-action interpolation tolerance
rank_tol_report  = [1e-4 1e-6 1e-8 1e-10]; % numerical-rank tolerances reported in section 4
r_max_common_basis = 80;   % largest common-basis rank tested in the projection-error sweep

n_validate_random = 10;    % final random held-out validation points (>=10 required)
rng_seed = 20260925;       % fixed seed for all random draws in this script

% Where cached exact snapshots and figures are stored. Both live under
% data/, which is already git-ignored (do not commit these).
cache_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_tikhonov_constant150_cache.mat');
fig_dir = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_tikhonov_figs');
if ~isfolder(fig_dir); mkdir(fig_dir); end
save_figs = true;

fprintf('R=%.3g, N_c=%d, N_f=%d, N_peanut=%d, n=2*N_c=%d\n', ...
    R,opt.N_c,opt.N_f,opt.N_peanut,2*opt.N_c);
fprintf('delta in [%.3e, %.3e] R, q_levels = %s\n', ...
    delta_min/R,delta_max/R,mat2str(q_levels));
fprintf('rank_trunc_tol=%.1e, matrix_action_tol=%.1e\n\n', ...
    rank_trunc_tol,matrix_action_tol);
fprintf('use_tikhonov=%d, tikhonov_tol=lambda/sigma_max=%.1e\n\n', ...
    opt.use_tikhonov,opt.tikhonov_tol);

%% ------------------------------------------------------------------
%  Fixed discretisation grids (independent of delta)
%  ------------------------------------------------------------------
N_c = opt.N_c; N_f = opt.N_f;
a_c = opt.a_c; a_f = opt.a_f;
Rp_c = opt.Rp_c; Rp_f = opt.Rp_f;

nout_c = ceil(a_c*N_c); tout_c = linspace(0,2*pi,nout_c+1)'; tout_c = tout_c(1:end-1);
grids.rbase_out_c = R*(cos(tout_c)+1i*sin(tout_c));
tin_c = linspace(0,2*pi,N_c+1)'; tin_c = tin_c(1:end-1);
grids.rbase_in_c = Rp_c*(cos(tin_c)+1i*sin(tin_c));
tin_f = linspace(0,2*pi,N_f+1)'; tin_f = tin_f(1:end-1);
grids.rbase_in_f = Rp_f*(cos(tin_f)+1i*sin(tin_f));
nout_f = ceil(a_f*N_f); tout_f = linspace(0,2*pi,nout_f+1)'; tout_f = tout_f(1:end-1);
grids.rout_base_f = R*(cos(tout_f)+1i*sin(tout_f));

%% ------------------------------------------------------------------
%  Self-test of the generic interpolation machinery (synthetic, fast)
%  ------------------------------------------------------------------
selfTestCapInterp();

%% ------------------------------------------------------------------
%  Exactly two comparison ranges, split at the former activation gap.
%  ------------------------------------------------------------------
delta_star = (R-Rp_f)^2/Rp_f;
assert(delta_min < delta_star && delta_star < delta_max, ...
    'The comparison split must lie strictly inside the interpolation range.');
fprintf('delta_star/R = %.8f (comparison split)\n',delta_star/R);

% With the ellipse discretisation frozen there is no branch at delta_star,
% so the two ranges share their endpoint exactly.
transition_margin = 0;
subintervals = struct( ...
    'name',{'small-gap range','large-gap range'}, ...
    'delta_lo',{delta_min,delta_star}, ...
    'delta_hi',{delta_star,delta_max}, ...
    'ellipse_active',{true,true});
fprintf(['Using two comparison ranges sharing delta_star/R=%.8f; ', ...
    'reference location=%s.\n\n'],delta_star/R,reference_location);

parametrisations = struct( ...
    'name',    {'s = log(delta/R)','alpha = acosh(1+delta/(2R))'}, ...
    'tag',     {'s','alpha'}, ...
    'to_param',{@(d) log(d/R), @(d) acosh(1+d/(2*R))}, ...
    'to_delta',{@(s) R*exp(s), @(a) 2*R*(cosh(a)-1)});

% On each (regime,parametrisation) interval, build exact matrices C(x_j)
% at Chebyshev-Lobatto nodes. The second-form barycentric interpolant is
%
%   C_hat(x) = sum_j [w_j/(x-x_j)] C(x_j) / sum_j [w_j/(x-x_j)].
%
% Thus every matrix element shares the same stable scalar weights; no
% elementwise fitting or matrix inversion occurs during interpolation.

%% ------------------------------------------------------------------
%  Full pipeline on each of the two subintervals x each parametrisation
%  ------------------------------------------------------------------
results = struct([]);
for si = 1:numel(subintervals)
    sub = subintervals(si);
    fprintf('=== Subinterval %d/%d: %s, delta/R in [%.3e, %.3e] ===\n', ...
        si,numel(subintervals),sub.name,sub.delta_lo/R,sub.delta_hi/R);
    for p = 1:numel(parametrisations)
        param = parametrisations(p);
        fprintf('--- parametrisation: %s ---\n',param.name);
        res = runGapInterpPipeline(sub,param,R,opt,grids,q_levels, ...
            rank_trunc_tol,matrix_action_tol,rank_tol_report,r_max_common_basis, ...
            n_validate_random,rng_seed,cache_file,reference_location);
        res.subinterval_name = sub.name;
        res.param_name = param.name;
        results = [results res]; %#ok<AGROW>
        fprintf('\n');
    end
end

%% ------------------------------------------------------------------
%  Plots
%  ------------------------------------------------------------------
plotGapInterpResults(results,q_levels,fig_dir,save_figs);

%% ------------------------------------------------------------------
%  Summary table
%  ------------------------------------------------------------------
fprintf('\n=== Summary table ===\n');
fprintf('%-28s %-14s %4s %6s %4s %6s %6s %10s %10s %10s %10s %10s\n', ...
    'subinterval','param','q','k(F)','m','r_L','r_R','fro_err','act_err', ...
    't_build','t_full_mv','t_fact_mv');
for i = 1:numel(results)
    r = results(i);
    fprintf('%-28s %-14s %4d %6d %4s %6d %6d %10.2e %10.2e %10.2e %10.2e %10.2e\n', ...
        r.subinterval_name,r.param_name,q_levels(end),r.snapshot_rank_k, ...
        'n/a',r.r_selected,r.r_selected,r.final_max_fro_err,r.final_max_action_err, ...
        r.cost.t_build,r.cost.t_matvec_full,r.cost.t_matvec_factorized);
end

fprintf(['\nColumns: q = finest training grid size; k(F) = numerical rank of the\n' ...
    'snapshot matrix F (parametric rank, section 6); m = AAA support-point\n' ...
    'count (n/a: no genuine QR-AAA implementation found, see below);\n' ...
    'r_L,r_R = common-basis rank used for the reduced core B(s); fro_err,\n' ...
    'act_err = max held-out Frobenius/spectral relative errors at the finest\n' ...
    'level and random test points; t_build = exact-construction time;\n' ...
    't_full_mv,t_fact_mv = time to apply the explicit vs factorised online\n' ...
    'representation to a vector.\n\n']);

%% ------------------------------------------------------------------
%  QR-AAA availability
%  ------------------------------------------------------------------
report_qraaa_availability();

%% ------------------------------------------------------------------
%  Final recommendation
%  ------------------------------------------------------------------
print_final_recommendation(results,delta_star,R,rank_trunc_tol,matrix_action_tol);

save(fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_tikhonov_results.mat'), ...
    'results','subintervals','parametrisations','delta_star', ...
    'transition_margin','q_levels','matrix_action_tol','rank_trunc_tol', ...
    'r_max_common_basis','reference_location','opt','R','-v7.3');
fprintf(['Results saved to data/', ...
    'sep25_capacitance_cmap_gap_interp_tikhonov_results.mat (git-ignored).\n']);


%% ==================================================================
%  Local helper functions (orchestration only; atomic numerical helpers
%  live in experiments/cap_interp/*.m and are reused here).
%  ==================================================================

function res = runGapInterpPipeline(sub,param,R,opt,grids,q_levels, ...
    rank_trunc_tol,matrix_action_tol,rank_tol_report,r_max_common_basis, ...
    n_validate_random,rng_seed,cache_file,reference_location)
%RUNGAPINTERPPIPELINE Sections 2-7 of the task, on one (subinterval,
% parametrisation) pair. See the header of the main script for the overall
% approach; this function is intentionally one place so all quantities in
% the summary table/plots for a given row are produced consistently.

param_lo = param.to_param(sub.delta_lo);
param_hi = param.to_param(sub.delta_hi);

%% Section 2-3: nested Chebyshev refinement with held-out validation
n_lev = numel(q_levels);
max_fro_err = nan(1,n_lev-1);
max_action_err = nan(1,n_lev-1);
effective_dof_fine_range = nan(1,n_lev);
effective_dof_peanut_range = nan(1,n_lev);
n_image_range = nan(2,n_lev);

nodes_final = [];
C_final = [];
info_final = [];

for i = 1:n_lev-1
    q_cur = q_levels(i);
    q_next = q_levels(i+1);
    nodes_cur = chebLobattoNodes(param_lo,param_hi,q_cur);
    nodes_next = chebLobattoNodes(param_lo,param_hi,q_next);
    held_out = new_nodes(nodes_next,nodes_cur,param_hi-param_lo);

    delta_cur = param.to_delta(nodes_cur);
    delta_held = param.to_delta(held_out);

    [C_cur,info_cur] = getOrBuildCmapSnapshots(delta_cur,R,opt,grids,cache_file);
    [C_held,~] = getOrBuildCmapSnapshots(delta_held,R,opt,grids,cache_file);

    w = chebBarycentricWeights(q_cur);
    Chat_held = evalMatrixChebBary(held_out,nodes_cur,w,C_cur);
    [fro_e,act_e] = held_out_errors(Chat_held,C_held);
    max_fro_err(i) = max(fro_e);
    max_action_err(i) = max(act_e);

    effective_dof_fine_range(i) = range_or_const([info_cur.effective_dof_fine]);
    effective_dof_peanut_range(i) = range_or_const([info_cur.effective_dof_peanut]);
    n_image_range(:,i) = [range_or_const(arrayfun(@(s) s.n_image(1),info_cur)); ...
                           range_or_const(arrayfun(@(s) s.n_image(2),info_cur))];

    fprintf('  q=%3d->%3d: %2d held-out pts, max fro err=%.3e, max action err=%.3e\n', ...
        q_cur,q_next,numel(held_out),max_fro_err(i),max_action_err(i));

    if i == n_lev-1
        nodes_final = nodes_next;
        [C_final,info_final] = getOrBuildCmapSnapshots(param.to_delta(nodes_next),R,opt,grids,cache_file);
    end
end
effective_dof_fine_range(end) = range_or_const([info_final.effective_dof_fine]);
effective_dof_peanut_range(end) = range_or_const([info_final.effective_dof_peanut]);
n_image_range(:,end) = [range_or_const(arrayfun(@(s) s.n_image(1),info_final)); ...
                         range_or_const(arrayfun(@(s) s.n_image(2),info_final))];

active_final = arrayfun(@(s) any(s.n_image>0),info_final);
assert(all(active_final==sub.ellipse_active), ...
    'Training nodes crossed the ellipse activation boundary.');
fprintf(['  Tikhonov effective-DOF spans: fine %.3f, peanut %.3f; ', ...
    'ellipse-node-count spans %d/%d.\n'], ...
    max(effective_dof_fine_range),max(effective_dof_peanut_range), ...
    max(n_image_range(1,:)),max(n_image_range(2,:)));
if any(n_image_range(:)>0)
    fprintf(['  NOTE: individual ellipse nodes enter/leave inside this regime; ', ...
        'they are intentionally validated through, not used as extra splits.\n']);
end

%% Final random-point validation (fixed seed, >=10 points)
rng(rng_seed);
param_rand = param_lo + (param_hi-param_lo)*rand(n_validate_random,1);
delta_rand = param.to_delta(param_rand);
[C_rand,~] = getOrBuildCmapSnapshots(delta_rand,R,opt,grids,cache_file);
w_final = chebBarycentricWeights(q_levels(end));
Chat_rand = evalMatrixChebBary(param_rand,nodes_final,w_final,C_final);
[fro_rand,act_rand] = held_out_errors(Chat_rand,C_rand);

res.q_levels = q_levels;
res.max_fro_err = max_fro_err;
res.max_action_err = max_action_err;
res.final_max_fro_err = max(max(fro_rand),max_fro_err(end));
res.final_max_action_err = max(max(act_rand),max_action_err(end));
res.random_param = param_rand;
res.random_fro_err = fro_rand;
res.random_action_err = act_rand;
res.effective_dof_fine_range = effective_dof_fine_range;
res.effective_dof_peanut_range = effective_dof_peanut_range;
res.n_image_range = n_image_range;
res.stopped_ok = res.final_max_action_err < matrix_action_tol;

fprintf('  Final random-point check (%d pts): max fro err=%.3e, max action err=%.3e -> %s\n', ...
    n_validate_random,max(fro_rand),max(act_rand), ...
    ternary(res.stopped_ok,'PASS (< matrix_action_tol)','FAIL (>= matrix_action_tol)'));

%% Section 4: individual + gap-dependent-part rank diagnostics at 3 gaps
param_3 = [param_lo; (param_lo+param_hi)/2; param_hi];
delta_3 = param.to_delta(param_3);
[C_3,~] = getOrBuildCmapSnapshots(delta_3,R,opt,grids,cache_file);
delta_ref = referenceDelta(sub.delta_lo,sub.delta_hi,reference_location);
[C_ref_stack,~] = getOrBuildCmapSnapshots( ...
    delta_ref,R,opt,grids,cache_file);
C_ref = C_ref_stack(:,:,1);
res.reference_location = reference_location;
res.reference_delta = delta_ref;

sigma_individual = zeros(size(C_3,1),3);
sigma_diff = zeros(size(C_3,1),3);
rank_individual = zeros(3,numel(rank_tol_report));
rank_diff = zeros(3,numel(rank_tol_report));
for k = 1:3
    sigma_individual(:,k) = svd(C_3(:,:,k));
    sigma_diff(:,k) = svd(C_3(:,:,k)-C_ref);
    rank_individual(k,:) = numericalRank(sigma_individual(:,k),rank_tol_report);
    rank_diff(k,:) = numericalRank(sigma_diff(:,k),rank_tol_report);
end
res.gap3_param = param_3;
res.gap3_delta = delta_3;
res.sigma_individual = sigma_individual;
res.sigma_diff = sigma_diff;
res.rank_individual = rank_individual;
res.rank_diff = rank_diff;
res.rank_tol_report = rank_tol_report;

fprintf(['  Reference delta/R=%.6g (%s).  Individual/Delta-A ranks at ', ...
    'rank_trunc_tol=%.1e: individual=[%d %d %d], Delta-A=[%d %d %d]\n'], ...
    delta_ref/R,reference_location, ...
    rank_trunc_tol,numericalRank(sigma_individual(:,1),rank_trunc_tol), ...
    numericalRank(sigma_individual(:,2),rank_trunc_tol),numericalRank(sigma_individual(:,3),rank_trunc_tol), ...
    numericalRank(sigma_diff(:,1),rank_trunc_tol),numericalRank(sigma_diff(:,2),rank_trunc_tol), ...
    numericalRank(sigma_diff(:,3),rank_trunc_tol));

%% Section 5: common spatial bases from the training snapshots
% "Held-out" here means genuinely unseen by the q=q_levels(end) training
% grid used to build U,V and B_train: only the final random points qualify.
dA_train = C_final - C_ref;
basis = commonBasisFromSnapshots(dA_train,r_max_common_basis);
dA_rand = C_rand - C_ref;

r_sweep = 1:r_max_common_basis;
proj_err_train = nan(size(r_sweep));
proj_err_held = nan(size(r_sweep));
for ir = 1:numel(r_sweep)
    r = r_sweep(ir);
    U_r = basis.U(:,1:r); V_r = basis.V(:,1:r);
    proj_err_train(ir) = max(proj_rel_err(dA_train,U_r,V_r,C_final));
    proj_err_held(ir) = max(proj_rel_err(dA_rand,U_r,V_r,C_rand));
end
r_selected = find(proj_err_held <= matrix_action_tol,1,'first');
spatial_reduction_pass = ~isempty(r_selected);
if ~spatial_reduction_pass
    r_selected = r_max_common_basis; % not fully converged; report as-is
end

res.basis_sigma_L = basis.sigma_L;
res.basis_sigma_R = basis.sigma_R;
res.r_sweep = r_sweep;
res.proj_err_train = proj_err_train;
res.proj_err_held = proj_err_held;
res.r_selected = r_selected;
res.spatial_reduction_pass = spatial_reduction_pass;

fprintf('  Common-basis rank r=%d gives max held-out projection error %.3e (tol %.1e): %s\n', ...
    r_selected,proj_err_held(min(r_selected,end)),matrix_action_tol, ...
    ternary(spatial_reduction_pass,'PASS','NOT CONVERGED BY r_max'));

%% Section 6: direct full-matrix vs reduced-core (B(s)) interpolation
U = basis.U(:,1:r_selected); V = basis.V(:,1:r_selected);
B_train = zeros(r_selected,r_selected,numel(nodes_final));
for j = 1:numel(nodes_final)
    B_train(:,:,j) = U'*dA_train(:,:,j)*V;
end
w_final = chebBarycentricWeights(q_levels(end));

% Evaluated only at the genuinely held-out random points, for the same
% reason as in section 5 above.
eval_pts = param_rand;
C_eval = C_rand;

Chat_direct = evalMatrixChebBary(eval_pts,nodes_final,w_final,C_final);
[fro_direct,act_direct] = held_out_errors(Chat_direct,C_eval);

Bhat = evalMatrixChebBary(eval_pts,nodes_final,w_final,B_train);
n_eval = numel(eval_pts);
Chat_reduced = zeros(size(C_ref,1),size(C_ref,2),n_eval);
for k = 1:n_eval
    Chat_reduced(:,:,k) = C_ref + U*Bhat(:,:,k)*V';
end
[fro_reduced,act_reduced] = held_out_errors(Chat_reduced,C_eval);

res.direct.max_fro_err = max(fro_direct);
res.direct.max_action_err = max(act_direct);
res.reduced.max_fro_err = max(fro_reduced);
res.reduced.max_action_err = max(act_reduced);

fprintf('  Direct full-matrix interpolation:  max fro=%.3e, max action=%.3e\n', ...
    res.direct.max_fro_err,res.direct.max_action_err);
fprintf('  Reduced-core (r=%d) interpolation: max fro=%.3e, max action=%.3e\n', ...
    r_selected,res.reduced.max_fro_err,res.reduced.max_action_err);

%% Snapshot (parametric) rank of F, section 6
n = size(dA_train,1);
q_train = numel(nodes_final);
F = reshape(dA_train,n*n,q_train).'; % q_train x n^2, F(j,:) = vec(dA(s_j))'
sigma_F = svd(F,'econ');
res.sigma_F = sigma_F;
res.snapshot_rank_k = numericalRank(sigma_F,rank_trunc_tol);
fprintf('  Snapshot matrix F: %d x %d, numerical rank k=%d at tol=%.1e\n', ...
    size(F,1),size(F,2),res.snapshot_rank_k,rank_trunc_tol);

%% Section 7: online cost comparison at one unseen gap
s_new = param_rand(1);
delta_new = param.to_delta(s_new);
n_rep = 200;

[~,info_new] = buildCanonicalPairCmap(delta_new,R,opt,grids); % exact build timing
t_build = info_new.t_build;

t0 = tic;
for r = 1:n_rep
    w_tmp = chebBarycentricWeights(q_levels(end)); %#ok<NASGU>
end
t_fit_weights = toc(t0)/n_rep; % (weights only; snapshots already cached)

g = randn(n,1);

t0 = tic;
for r = 1:n_rep
    Chat_new = evalMatrixChebBary(s_new,nodes_final,w_final,C_final); %#ok<NASGU>
end
t_eval_full = toc(t0)/n_rep;
Chat_new = evalMatrixChebBary(s_new,nodes_final,w_final,C_final);

t0 = tic;
for r = 1:n_rep
    y_full = Chat_new*g; %#ok<NASGU>
end
t_matvec_full = toc(t0)/n_rep;

t0 = tic;
for r = 1:n_rep
    Bhat_new = evalMatrixChebBary(s_new,nodes_final,w_final,B_train); %#ok<NASGU>
end
t_eval_reduced = toc(t0)/n_rep;
Bhat_new = evalMatrixChebBary(s_new,nodes_final,w_final,B_train);

t0 = tic;
for r = 1:n_rep
    y_fact = C_ref*g + U*(Bhat_new*(V'*g)); %#ok<NASGU>
end
t_matvec_factorized = toc(t0)/n_rep;

bytes_full_storage = q_train*n*n*8;
% The reduced representation must retain the dense reference C_ref. Count
% it in both storage and apply cost; omitting it would make the reduction
% look artificially cheaper for a standalone Cmap application.
bytes_factorized_storage = (n*n + n*r_selected*2 + ...
    q_train*r_selected*r_selected)*8;

res.cost.t_build = t_build;
res.cost.t_fit_weights = t_fit_weights;
res.cost.t_eval_full = t_eval_full;
res.cost.t_eval_reduced = t_eval_reduced;
res.cost.t_matvec_full = t_matvec_full;
res.cost.t_matvec_factorized = t_matvec_factorized;
res.cost.bytes_full_storage = bytes_full_storage;
res.cost.bytes_factorized_storage = bytes_factorized_storage;
res.cost.n = n;
res.cost.r = r_selected;
res.cost.dense_flop_estimate = n^2;
res.cost.increment_flop_estimate = 2*n*r_selected + r_selected^2;
res.cost.factorized_flop_estimate = n^2 + res.cost.increment_flop_estimate;

fprintf(['  Online cost @ unseen gap: build=%.2es, eval(full)=%.2es, mv(full)=%.2es, ' ...
    'mv(factorized,r=%d)=%.2es\n'],t_build,t_eval_full,t_matvec_full,r_selected,t_matvec_factorized);
fprintf(['  Apply cost: dense n^2=%d vs C_ref+reduced increment ', ...
    'n^2+2nr+r^2=%d (dense/factorized %.2gx); storage ratio %.2gx.\n\n'], ...
    res.cost.dense_flop_estimate,res.cost.factorized_flop_estimate, ...
    res.cost.dense_flop_estimate/res.cost.factorized_flop_estimate, ...
    bytes_full_storage/bytes_factorized_storage);

end

function held_out = new_nodes(nodes_next,nodes_cur,scale)
d = min(abs(nodes_next(:)-nodes_cur(:)'),[],2);
held_out = nodes_next(d > 1e-9*max(1,abs(scale)));
end

function v = range_or_const(x)
v = max(x)-min(x);
end

function s = ternary(cond,a,b)
if cond, s = a; else, s = b; end
end

function [fro_e,act_e] = held_out_errors(Chat,C)
% Chat,C: n x n x m stacks. fro_e(k) = ||C_k-Chat_k||_F/||C_k||_F,
% act_e(k) = ||C_k-Chat_k||_2/||C_k||_2 (exact matrix 2-norm, i.e. exactly
% the max-over-unit-vectors quantity in the task description).
m = size(C,3);
fro_e = zeros(m,1); act_e = zeros(m,1);
for k = 1:m
    D = Chat(:,:,k)-C(:,:,k);
    fro_e(k) = norm(D,'fro')/max(norm(C(:,:,k),'fro'),eps);
    act_e(k) = norm(D,2)/max(norm(C(:,:,k),2),eps);
end
end

function e = proj_rel_err(dA,U_r,V_r,C_for_norm)
% e(k) = ||dA_k - U_r U_r' dA_k V_r V_r'||_2 / ||C_for_norm_k||_2
m = size(dA,3);
e = zeros(m,1);
for k = 1:m
    P = U_r*(U_r'*dA(:,:,k)*V_r)*V_r';
    e(k) = norm(dA(:,:,k)-P,2)/max(norm(C_for_norm(:,:,k),2),eps);
end
end

function report_qraaa_availability()
fprintf('=== QR-AAA availability ===\n');
has_setvalued = exist('qraaa','file')==2 || exist('setAAA','file')==2 || ...
    exist('vectorAAA','file')==2 || exist('matrixAAA','file')==2;
has_scalar_aaa = exist('aaa','file')==2;
if has_setvalued
    fprintf(['A matrix/set-valued AAA implementation was found on the MATLAB\n' ...
        'path and should be used for the QR-AAA comparison (not implemented\n' ...
        'automatically here; wire it in at the marked location).\n\n']);
    return
end

fprintf(['No matrix- or set-valued (QR-AAA) implementation was found on\n' ...
    'the repository or MATLAB path.\n']);
if has_scalar_aaa
    fprintf(['A scalar, entrywise aaa.m IS on the path, but it only fits one\n' ...
        'function at a time; running it once per matrix entry is not QR-AAA\n' ...
        '(it ignores the shared pole/support structure across entries, and\n' ...
        'its cost scales with n^2 independent scalar fits). It is therefore\n' ...
        'NOT substituted here.\n']);
else
    fprintf('No scalar aaa.m was found either.\n');
end
fprintf(['To add QR-AAA correctly this script already provides the required\n' ...
    'interface: F (q x n^2, snapshot rows = vec(Delta A(s_j))^T) and its\n' ...
    'economical SVD/QR (snapshot rank k, already computed above as\n' ...
    'res.snapshot_rank_k). The missing piece is a genuine set-valued AAA\n' ...
    'solver that takes the k dominant row combinations of F (a length-q set\n' ...
    'of k scalar functions of s) and returns one shared barycentric\n' ...
    'support-point set for all of them; the reconstructed matrix interpolant\n' ...
    'would then reuse the QR factor to map back to n^2 entries. No such\n' ...
    'solver exists in this repository or on this MATLAB path at the time of\n' ...
    'writing, so it was not fitted; the Chebyshev and common-basis\n' ...
    '(reduced-core) results above stand in as the completed, reliable\n' ...
    'baseline and reduced-order methods.\n\n']);
end

function plotGapInterpResults(results,q_levels,fig_dir,save_figs)
%PLOTGAPINTERPRESULTS All required plots for one script run: individual and
% variation (Delta A) singular values, common-basis singular values,
% snapshot singular values, error vs q, held-out error vs s/alpha,
% projection error vs rank. One figure per (subinterval,parametrisation)
% result.

for i = 1:numel(results)
    r = results(i);
    tag = sprintf('%s | %s',r.subinterval_name,r.param_name);

    f = figure('Name',['Individual and variation singular values: ' tag],'Visible','on');
    subplot(1,2,1);
    semilogy(r.sigma_individual./r.sigma_individual(1,:),'-o','MarkerSize',3);
    xlabel('j'); ylabel('\sigma_j(A)/\sigma_1(A)');
    legend({'smallest gap','coordinate midpoint','largest gap'},'Location','best');
    title(['A(s): ' tag]); grid on;
    subplot(1,2,2);
    semilogy(r.sigma_diff./max(r.sigma_diff(1,:),eps),'-o','MarkerSize',3);
    xlabel('j'); ylabel('\sigma_j(\Delta A)/\sigma_1(\Delta A)');
    legend({'smallest gap','mid gap','largest gap'},'Location','best');
    title(['\Delta A(s) = A(s)-A_{ref}: ' tag]); grid on;
    save_fig(f,fig_dir,sprintf('sigma_individual_diff_%02d',i),save_figs);

    f = figure('Name',['Common-basis and snapshot singular values: ' tag],'Visible','on');
    subplot(1,2,1);
    semilogy(r.basis_sigma_L,'-o','DisplayName','left (U)','MarkerSize',3); hold on;
    semilogy(r.basis_sigma_R,'-s','DisplayName','right (V)','MarkerSize',3);
    xlabel('r'); ylabel('\sigma_r'); legend('Location','best');
    title(['Common-basis singular values: ' tag]); grid on;
    subplot(1,2,2);
    semilogy(r.sigma_F,'-o','MarkerSize',3);
    xlabel('j'); ylabel('\sigma_j(F)');
    title(['Snapshot matrix F singular values: ' tag]); grid on;
    save_fig(f,fig_dir,sprintf('sigma_basis_snapshot_%02d',i),save_figs);

    f = figure('Name',['Projection error vs rank: ' tag],'Visible','on');
    semilogy(r.r_sweep,r.proj_err_train,'-o','DisplayName','training gaps','MarkerSize',3); hold on;
    semilogy(r.r_sweep,r.proj_err_held,'-s','DisplayName','held-out gaps','MarkerSize',3);
    xline(r.r_selected,'--','DisplayName',sprintf('selected r=%d',r.r_selected));
    xlabel('common-basis rank r'); ylabel('max_g ||\Delta A - UU^*\Delta A VV^*||_2/||A||_2');
    legend('Location','best'); title(['Projection error vs rank: ' tag]); grid on;
    save_fig(f,fig_dir,sprintf('proj_err_vs_rank_%02d',i),save_figs);

    f = figure('Name',['Interpolation error vs q: ' tag],'Visible','on');
    subplot(1,2,1);
    semilogy(q_levels(2:end),r.max_fro_err,'-o','DisplayName','Frobenius'); hold on;
    semilogy(q_levels(2:end),r.max_action_err,'-s','DisplayName','spectral/action');
    xlabel('q (finer level)'); ylabel('max held-out relative error');
    legend('Location','best'); title(['Held-out error vs q: ' tag]); grid on;
    subplot(1,2,2);
    semilogy(r.random_param,r.random_fro_err,'o','DisplayName','Frobenius'); hold on;
    semilogy(r.random_param,r.random_action_err,'s','DisplayName','spectral/action');
    xlabel('gap parameter (random test points; s or \alpha, may be negative)'); ylabel('relative error');
    legend('Location','best'); title(['Final random-point errors: ' tag]); grid on;
    save_fig(f,fig_dir,sprintf('err_vs_q_and_s_%02d',i),save_figs);
end

end

function save_fig(f,fig_dir,name,save_figs)
% Keep every figure visible; save_figs controls only whether a copy is saved.
if save_figs
    saveas(f,fullfile(fig_dir,[name '.png']));
end
drawnow;
end

function print_final_recommendation(results,delta_star,R,rank_trunc_tol,matrix_action_tol)
fprintf('\n=== Final recommendation ===\n');

all_pass = all(arrayfun(@(r) r.stopped_ok,results));
max_act = max(arrayfun(@(r) r.final_max_action_err,results));
fprintf(['1. Barycentric interpolation across the two comparison ranges: %s. ', ...
    'The largest held-out action error is %.2e (target %.1e). No ', ...
    'rank-based subintervals are used.\n'], ...
    ternary_str(all_pass,'PASS','NOT FULLY CONVERGED'),max_act,matrix_action_tol);

q_needed = arrayfun(@(r) first_q_below_tol(r,matrix_action_tol),results);
q_finite = q_needed(isfinite(q_needed));
if isempty(q_finite)
    fprintf(['2. No tested grid through q=%d reached the %.1e action-error ', ...
        'target on every subinterval; see err_vs_q_and_s_*.png.\n'], ...
        results(1).q_levels(end),matrix_action_tol);
else
    q_max = max(q_finite);
    fprintf(['2. Up to q=%.0f exact Cmap snapshots were needed on each ', ...
        'successful subinterval. See err_vs_q_and_s_*.png.\n'],q_max);
end

r_report = arrayfun(@(r) r.r_selected,results);
k_report = arrayfun(@(r) r.snapshot_rank_k,results);
spatial_pass = arrayfun(@(r) r.spatial_reduction_pass,results);
storage_ratio = arrayfun(@(r) ...
    r.cost.bytes_full_storage/r.cost.bytes_factorized_storage,results);
apply_ratio = arrayfun(@(r) ...
    r.cost.t_matvec_full/r.cost.t_matvec_factorized,results);
fprintf(['3. Spatial reduction of Delta C=C-C_ref at rank tolerance %.1e: ', ...
    'r=%s, snapshot ranks k=%s, and the accuracy target is met in %d/%d ', ...
    'cases. Including dense C_ref, storage is reduced by %.2g--%.2gx; ', ...
    'measured dense/factorized matvec ratios are %.2g--%.2gx (>1 is ', ...
    'faster). It is therefore primarily a parametric-storage reduction ', ...
    'unless C_ref is already available or applied separately.\n'], ...
    rank_trunc_tol,mat2str(r_report),mat2str(k_report),sum(spatial_pass), ...
    numel(spatial_pass),min(storage_ratio),max(storage_ratio), ...
    min(apply_ratio),max(apply_ratio));

fprintf(['4. Interval policy: one split at delta_star/R=%.8f gives exactly ', ...
    'two comparison ranges. Tikhonov removes hard pseudoinverse-rank ', ...
    'switches, and ellipse_constant=true freezes the 150-node candidate ', ...
    'discretisation throughout both ranges.\n'],delta_star/R);

s_err = arrayfun(@(r) r.final_max_action_err, ...
    results(strcmp({results.param_name},'s = log(delta/R)')));
a_err = arrayfun(@(r) r.final_max_action_err, ...
    results(strcmp({results.param_name},'alpha = acosh(1+delta/(2R))')));
mean_s = mean(s_err);
mean_a = mean(a_err);
relative_difference = abs(mean_a-mean_s)/min(mean_a,mean_s);
if relative_difference < 0.1
    better = 'no material difference between the two coordinates';
elseif mean_a < mean_s
    better = 'alpha = acosh(1+delta/(2R))';
else
    better = 's = log(delta/R)';
end
fprintf(['5. Mean final action error: log(delta/R) %.2e, alpha %.2e. ', ...
    'Recommendation from this run: %s.\n'],mean_s,mean_a,better);

end

function s = ternary_str(cond,a,b)
if cond, s = a; else, s = b; end
end

function q = first_q_below_tol(r,tol)
i = find(r.max_action_err < tol,1,'first');
if isempty(i)
    q = NaN;
else
    q = r.q_levels(i+1);
end
end

function delta_ref = referenceDelta(delta_lo,delta_hi,location)
switch location
    case 'minimum'
        delta_ref = delta_lo;
    case 'midpoint'
        delta_ref = (delta_lo+delta_hi)/2;
    case 'maximum'
        delta_ref = delta_hi;
    otherwise
        error('Unknown reference location "%s".',location);
end
end
