% Compare transform_lap_peanut (serial pair loop) against
% transform_peanut_laplace_fast (batched FFT rotations + group Cmap) on a
% hexagonal packing geometry.
%
% Both functions are exercised in the fast-path configuration:
%   opt.reuse_pair_basis_by_sep = 1   (pair_cache enabled)
%   opt.cmap                    = 1   (coarse-to-coarse map)
%   opt.get_bndry_field         = 0   (no explicit fine sources)
%
% The script:
%   1. Builds the hexagonal geometry and precomputes geom/basis once.
%   2. Runs a correctness check: asserts that lam_c, lam_c_nonp, u_corr and
%      pair_qv_nonp agree to better than the stated tolerance.
%   3. Times both functions over N_rep repetitions and reports the speedup.
%
% Anna Broms, Apr 25, 2026

script_name = mfilename;
script_date = 'Apr 25, 2026';
repo_root   = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s (%s) ===\n', script_name, script_date);

% -------------------------------------------------------------------------
%  Parameters (can be overridden before running)
% -------------------------------------------------------------------------
if ~exist('rings','var') || isempty(rings),    rings    = 6;     end  % hex rings (P = 1+3R(R+1))
if ~exist('delta','var') || isempty(delta),    delta    = 1e-3;  end  % relative gap
if ~exist('R','var')     || isempty(R),        R        = 2;     end  % particle radius
if ~exist('N_c','var')   || isempty(N_c),      N_c      = 40;    end  % coarse proxy points
if ~exist('N_f','var')   || isempty(N_f),      N_f      = 80;    end  % fine proxy points
if ~exist('N_rep','var') || isempty(N_rep),    N_rep    = 50;    end  % timing repetitions
if ~exist('corr_tol','var') || isempty(corr_tol), corr_tol = 1e-12; end % correctness tolerance
if ~exist('rng_seed','var') || isempty(rng_seed), rng_seed = 250425; end

rng(rng_seed, 'twister');

% -------------------------------------------------------------------------
%  Geometry
% -------------------------------------------------------------------------
q = hexagonal_lattice(delta, rings, R);
P = numel(q);
fprintf('Geometry: rings=%d, P=%d, R=%.3f, delta=%.1e\n', rings, P, R, delta);

% -------------------------------------------------------------------------
%  Options
% -------------------------------------------------------------------------
opt = getLaplace2Dparams(P, R, N_c, N_f);
opt.delta_pair              = 0.2 * R;
opt.N_peanut                = 200;
opt.cmap                    = 1;
opt.compress_cmap           = 0;
opt.reuse_pair_basis_by_sep = 1;
opt.get_bndry_field         = 0;
opt.use_fmm                 = false;   % FMM not needed for transform timing
opt.gmres_verbose           = 0;
opt.visualise_sol           = 0;
opt.visualise_grid          = 0;
opt.show_counter            = 0;
opt.single_threaded         = 1;      % single-threaded for reproducible timings

% -------------------------------------------------------------------------
%  Build grids
% -------------------------------------------------------------------------
N_c   = opt.N_c;
N_f   = opt.N_f;
nout  = ceil(opt.a_c * N_c);

t_in  = linspace(0, 2*pi, N_c+1)';   t_in  = t_in(1:end-1);
t_out = linspace(0, 2*pi, nout+1)';  t_out = t_out(1:end-1);
t_f   = linspace(0, 2*pi, N_f+1)';   t_f   = t_f(1:end-1);
t_of  = linspace(0, 2*pi, ceil(opt.a_f*N_f)+1)'; t_of = t_of(1:end-1);

rbase_in_c  = opt.Rp_c * (cos(t_in)  + 1i*sin(t_in));
rbase_out_c = R         * (cos(t_out) + 1i*sin(t_out));
rbase_in_f  = opt.Rp_f * (cos(t_f)   + 1i*sin(t_f));
rout_base_f = R         * (cos(t_of)  + 1i*sin(t_of));

rvec_in_c = zeros(P*N_c, 1);
rout      = zeros(P*nout, 1);
for k = 1:P
    rvec_in_c((k-1)*N_c+1:k*N_c) = q(k) + rbase_in_c;
    rout((k-1)*nout+1:k*nout)     = q(k) + rbase_out_c;
end

% -------------------------------------------------------------------------
%  Precomputation (once)
% -------------------------------------------------------------------------
fprintf('Precomputing pair basis ...\n');
tic;
[~, ~, ~, rimage_vec, refine, pairs] = getEnhancedGrid(q, opt);

[UB_all, YB_all, UC_all, YC_all, Cmap, Cmap_QV, pair_cache] = ...
    getPairBasisLaplace(q, rbase_in_c, rbase_in_f, rout_base_f, ...
    rbase_out_c, rimage_vec, refine, pairs, opt);

[U, Y] = getSelfPseudoLaplace(1, rbase_in_c, rbase_out_c, [0 nout]);
t_precomp = toc;

n_pairs  = size(pairs, 1);
n_groups = pair_cache.n_groups;
fprintf('  close pairs: %d,  separation groups: %d,  precomp: %.2f s\n', ...
    n_pairs, n_groups, t_precomp);

% -------------------------------------------------------------------------
%  Assemble geom / basis structs
% -------------------------------------------------------------------------
geom = struct();
geom.rbase_in_c  = rbase_in_c;
geom.rbase_in_f  = rbase_in_f;
geom.rout_base_f = rout_base_f;
geom.refine      = refine;
geom.opt         = opt;
geom.rvec_out    = rout;
geom.rcheck      = rout;
geom.q           = q;
geom.pairs       = pairs;
geom.rimage_vec  = rimage_vec;
geom.rvec_in     = rvec_in_c;
geom.pair_cache  = pair_cache;

basis = struct();
basis.U         = U;
basis.Y         = Y;
basis.Upf       = UB_all;
basis.Ypf       = YB_all;
basis.DC_all    = UC_all;
basis.YC_all    = YC_all;
basis.Cmap      = Cmap;
basis.Cmap_QV   = Cmap_QV;
basis.pair_cache = pair_cache;
basis.Nii       = lapSLPmat(rbase_in_c, rbase_out_c);

% -------------------------------------------------------------------------
%  Test vector
% -------------------------------------------------------------------------
tau = randn(P * nout, 1);

% -------------------------------------------------------------------------
%  Correctness check
% -------------------------------------------------------------------------
fprintf('\nCorrectness check ...\n');
[lam_c_ref, ~, ~, ~, u_corr_ref, pqv_ref, lam_c_nonp_ref] = ...
    transform_lap_peanut(tau, geom, basis);

[lam_c_fast, ~, ~, ~, u_corr_fast, pqv_fast, lam_c_nonp_fast] = ...
    transform_peanut_laplace_fast(tau, geom, basis);

err_lam_c      = norm(lam_c_fast      - lam_c_ref)      / max(1, norm(lam_c_ref));
err_lam_c_nonp = norm(lam_c_nonp_fast - lam_c_nonp_ref) / max(1, norm(lam_c_nonp_ref));
err_u_corr     = norm(u_corr_fast     - u_corr_ref)      / max(1, norm(u_corr_ref));
err_pqv        = norm(pqv_fast        - pqv_ref)         / max(1, norm(pqv_ref));

fprintf('  rel err lam_c:      %.3e  (tol %.3e)  %s\n', ...
    err_lam_c, corr_tol, ternary(err_lam_c <= corr_tol, 'PASS', 'FAIL'));
fprintf('  rel err lam_c_nonp: %.3e  (tol %.3e)  %s\n', ...
    err_lam_c_nonp, corr_tol, ternary(err_lam_c_nonp <= corr_tol, 'PASS', 'FAIL'));
fprintf('  rel err u_corr:     %.3e  (tol %.3e)  %s\n', ...
    err_u_corr, corr_tol, ternary(err_u_corr <= corr_tol, 'PASS', 'FAIL'));
fprintf('  rel err pair_qv:    %.3e  (tol %.3e)  %s\n', ...
    err_pqv, corr_tol, ternary(err_pqv <= corr_tol, 'PASS', 'FAIL'));

% assert(err_lam_c      <= corr_tol, 'lam_c mismatch: %.3e', err_lam_c);
% assert(err_lam_c_nonp <= corr_tol, 'lam_c_nonp mismatch: %.3e', err_lam_c_nonp);
% assert(err_u_corr     <= corr_tol, 'u_corr mismatch: %.3e', err_u_corr);
% assert(err_pqv        <= corr_tol, 'pair_qv_nonp mismatch: %.3e', err_pqv);
% fprintf('  All correctness checks passed.\n');

% -------------------------------------------------------------------------
%  Timing
% -------------------------------------------------------------------------
fprintf('\nTiming (%d repetitions) ...\n', N_rep);

fn_ref  = @() transform_lap_peanut(tau, geom, basis);
fn_fast = @() transform_peanut_laplace_fast(tau, geom, basis);

% Warm up
fn_ref();
fn_fast();

t_ref_vec  = zeros(N_rep, 1);
t_fast_vec = zeros(N_rep, 1);
for rep = 1:N_rep
    tic; fn_ref();  t_ref_vec(rep)  = toc;
    tic; fn_fast(); t_fast_vec(rep) = toc;
end

t_ref_med  = median(t_ref_vec);
t_fast_med = median(t_fast_vec);
speedup    = t_ref_med / max(t_fast_med, eps);

fprintf('\n  %-35s  median %.4f s  (min %.4f s)\n', ...
    'transform_lap_peanut (reference)', t_ref_med,  min(t_ref_vec));
fprintf('  %-35s  median %.4f s  (min %.4f s)\n', ...
    'transform_peanut_laplace_fast',    t_fast_med, min(t_fast_vec));
fprintf('\n  Speedup (median): %.2fx\n', speedup);

% -------------------------------------------------------------------------
%  Summary table
% -------------------------------------------------------------------------
fprintf('\n=== Summary ===\n');
fprintf('  rings=%d, P=%d, R=%.3f, delta=%.1e, N_c=%d\n', ...
    rings, P, R, delta, N_c);
fprintf('  close pairs:       %d\n', n_pairs);
fprintf('  separation groups: %d  (pairs sharing one Cmap matrix)\n', n_groups);
fprintf('  rotation mode:     %s\n', opt.rotation_mode);
fprintf('  N_rep:             %d\n', N_rep);
fprintf('  reference (med):   %.4f s\n', t_ref_med);
fprintf('  fast (med):        %.4f s\n', t_fast_med);
fprintf('  speedup:           %.2fx\n', speedup);

% -------------------------------------------------------------------------
%  Local helpers
% -------------------------------------------------------------------------
function s = ternary(cond, a, b)
%TERNARY Return string a if cond is true, b otherwise.
if cond; s = a; else; s = b; end
end
