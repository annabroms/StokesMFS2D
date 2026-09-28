function [lam_c,lam_self,lam_f,lam_e,u_corr,pair_qv_nonp,...
        lam_c_nonp,lam_self_nonp,lam_f_nonp,lam_e_nonp] = ...
    transform_peanut_laplace_fast(tau, geom, basis)
%TRANSFORM_PEANUT_LAPLACE_FAST Fast variant of TRANSFORM_LAP_PEANUT for the
%  case opt.reuse_pair_basis_by_sep=1 and opt.cmap=1.
%
% Syntax:
%   [lam_c,lam_self,lam_f,lam_e,u_corr] = ...
%       transform_peanut_laplace_fast(tau,geom,basis)
%   [lam_c,lam_self,lam_f,lam_e,u_corr,pair_qv_nonp,lam_c_nonp, ...
%       lam_self_nonp,lam_f_nonp,lam_e_nonp] = ...
%       transform_peanut_laplace_fast(tau,geom,basis)
%
% The function is a drop-in replacement for TRANSFORM_LAP_PEANUT restricted
% to the fast path (pair_cache.enabled=true, cmap=true,
% get_bndry_field=false).  Fine-source outputs (lam_f, lam_e, lam_f_nonp,
% lam_e_nonp) are always returned empty because explicit pair sources are
% not required in this path.
%
% Speedups over TRANSFORM_LAP_PEANUT (pair loop):
%   1. Batch FFT phase-rotations: all pair inputs / outputs are gathered
%      into (N_c x 2*n_pairs) matrices and rotated with a single FFT call.
%      Supports numeric phases (fft mode), struct specs with mode='fft' or
%      'shift', and mode='oversampled_fft' via batch zero-padded FFT.
%   2. Group-batched Cmap: pairs sharing the same canonical separation are
%      processed as one matrix-matrix multiply per group instead of n_pairs
%      separate mat-vecs.  On regular geometries (e.g., hexagonal packing)
%      most pairs share a single group, giving near-optimal reuse.
%   3. Vectorised scatter: pair corrections are accumulated with accumarray
%      instead of per-pair indexed assignment.
%   4. u_corr: computed in a serial loop (per-pair geometry differs);
%      per-pair data is read from the batch arrays rather than re-extracted.
%
% Inputs/outputs identical to TRANSFORM_LAP_PEANUT.
%
% See also: transform_lap_peanut, getPairBasisLaplace, getPeanutBlockLaplace,
%   matvec_lap_peanut_enhanced.
%
% Anna Broms, Apr 2026

%% Validate fast-path conditions
pair_cache = basis.pair_cache;
opt        = geom.opt;

if ~strcmp(getLaplaceInterpolationMode(opt),'none')
    error('transform_peanut_laplace_fast:InterpolationUnsupported', ...
        ['Interpolated maps use the general transform so reduced cores ', ...
         'remain factored.']);
end

if ~pair_cache.enabled
    error('transform_peanut_laplace_fast:NoPairCache', ...
        ['opt.reuse_pair_basis_by_sep must be true (pair_cache.enabled=1) ', ...
         'for the fast path. Use transform_lap_peanut for the general case.']);
end

use_cmap = isfield(opt,'cmap') && opt.cmap;
if ~use_cmap
    error('transform_peanut_laplace_fast:NoCmap', ...
        ['opt.cmap must be true for the fast path. ', ...
         'Use transform_lap_peanut for the general case.']);
end

get_bndry_field = logical(getOptField(opt,'get_bndry_field',true));
if get_bndry_field
    error('transform_peanut_laplace_fast:BndryFieldUnsupported', ...
        ['opt.get_bndry_field=true requires explicit fine sources, which are ', ...
         'not computed in the fast path. Set opt.get_bndry_field=0 or use ', ...
         'transform_lap_peanut instead.']);
end

%% Unpack geometry
q          = geom.q;
pairs      = geom.pairs;
rvec_out   = geom.rvec_out;
rbase_in_c = geom.rbase_in_c;
if isfield(geom,'rcheck') && ~isempty(geom.rcheck)
    rcheck_out = geom.rcheck;
else
    rcheck_out = rvec_out;
end

P       = numel(q);
N_c     = opt.N_c;
N_large = numel(rvec_out) / P;
N_check = numel(rcheck_out) / P;

if isfield(opt,'project_charge') && ~isempty(opt.project_charge)
    project_charge = logical(opt.project_charge);
else
    project_charge = false;
end

%% Fine-source outputs (never needed in this path)
lam_f      = [];
lam_e      = [];
lam_f_nonp = [];
lam_e_nonp = [];

%% Phase 1: one-body coarse mapping
% All one-body blocks are identical, so batch them as particle columns
% while keeping the U{1} then Y{1} ordering for a backward-stable apply.
tau_blocks      = reshape(tau(1:P*N_large), N_large, P);
lam_blocks_nonp = basis.Y{1} * (basis.U{1} * tau_blocks);  % (N_c x P)
if project_charge
    lam_blocks = lam_blocks_nonp - mean(lam_blocks_nonp, 1);
else
    lam_blocks = lam_blocks_nonp;
end
lam_c       = reshape(lam_blocks,     [], 1);
lam_c_nonp  = reshape(lam_blocks_nonp, [], 1);
lam_self      = lam_c;
lam_self_nonp = lam_c_nonp;

pair_qv_nonp = zeros(P, 1);
u_corr       = zeros(P * N_check, 1);

n_pairs = size(pairs, 1);
if n_pairs == 0
    return
end

%% Phase 2: batch pair processing

% --- 2a. Vectorised gather: one (N_c x n_pairs) block per particle -------
%  idx_i_mat(k,r) = global coarse-source index for mode k of particle
%  pairs(r,1), and similarly for idx_j_mat and pairs(r,2).
local     = (1:N_c)';
idx_i_mat = bsxfun(@plus, (pairs(:,1)-1)' * N_c, local);  % N_c x n_pairs
idx_j_mat = bsxfun(@plus, (pairs(:,2)-1)' * N_c, local);  % N_c x n_pairs

pair_top = lam_self(idx_i_mat);   % N_c x n_pairs  (particle 1)
pair_bot = lam_self(idx_j_mat);   % N_c x n_pairs  (particle 2)

% --- 2b. Detect non-identity rotations ------------------------------------
metas     = pair_cache.meta;
rots      = [metas.rot];
needs_rot = abs(rots - 1) > 100 * eps(max(1, abs(rots)));
rot_rows  = find(needs_rot);

% --- 2c. Batch forward rotation to canonical frame -----------------------
[pair_top, pair_bot] = batchRotatePair(pair_top, pair_bot, metas, ...
    rot_rows, N_c, isreal(tau), 'phase_c');

% --- 2d. Group-batched Cmap + Cmap_QV ------------------------------------
%  pair_rhs_mat(:,r) = [rotated lam_i; rotated lam_p2] for pair r.
pair_rhs_mat = [pair_top; pair_bot];     % (2*N_c x n_pairs)
tau_nonp_mat = zeros(2*N_c, n_pairs);
qv_mat       = zeros(2, n_pairs);

group_ids   = [metas.group_id];
unique_gids = unique(group_ids);
for k = 1:numel(unique_gids)
    gid  = unique_gids(k);
    cols = find(group_ids == gid);
    g    = pair_cache.groups(gid);
    % Matrix-matrix multiply: one Cmap apply for all pairs in this group.
    tau_nonp_mat(:, cols) = g.Cmap    * pair_rhs_mat(:, cols);
    qv_mat(:,       cols) = g.Cmap_QV * pair_rhs_mat(:, cols);
end

% --- 2e. Batch inverse rotation to lab frame -----------------------------
tau_nonp_top = tau_nonp_mat(1:N_c,     :);   % N_c x n_pairs  (particle 1)
tau_nonp_bot = tau_nonp_mat(N_c+1:end, :);   % N_c x n_pairs  (particle 2)
[tau_nonp_top, tau_nonp_bot] = batchRotatePair(tau_nonp_top, tau_nonp_bot, ...
    metas, rot_rows, N_c, isreal(tau), 'phase_c_inv');

% --- 2f. Charge projection (column-wise mean subtraction) ----------------
if project_charge
    tau_proj_top = tau_nonp_top - mean(tau_nonp_top, 1);
    tau_proj_bot = tau_nonp_bot - mean(tau_nonp_bot, 1);
else
    tau_proj_top = tau_nonp_top;
    tau_proj_bot = tau_nonp_bot;
end

% --- 2g. Vectorised scatter -----------------------------------------------
%  Flatten scatter-row indices and values, then accumulate in one call.
scatter_rows = [idx_i_mat(:); idx_j_mat(:)];          % (2*N_c*n_pairs x 1)
n_coarse     = P * N_c;

delta_c      = accumarray(scatter_rows, ...
    [tau_proj_top(:);  tau_proj_bot(:)],  [n_coarse, 1]);
delta_c_nonp = accumarray(scatter_rows, ...
    [tau_nonp_top(:); tau_nonp_bot(:)], [n_coarse, 1]);

lam_c      = lam_self      + delta_c;
lam_c_nonp = lam_self_nonp + delta_c_nonp;

qv_rows      = [pairs(:,1); pairs(:,2)];               % (2*n_pairs x 1)
pair_qv_nonp = accumarray(qv_rows, qv_mat(:), [P, 1]);

%% Phase 3: u_corr (serial loop; per-pair data read from batch arrays)
use_ucorr_colloc = isequal(rcheck_out, rvec_out);

for row = 1:n_pairs
    i   = pairs(row, 1);
    p2  = pairs(row, 2);
    meta = metas(row);

    block_i  = (i-1)*N_check+1 : i*N_check;
    block_p2 = (p2-1)*N_check+1 : p2*N_check;
    rout_pair  = [rcheck_out(block_i); rcheck_out(block_p2)];
    rin_pair_c = [q(i)+rbase_in_c; q(p2)+rbase_in_c];

    % Per-pair quantities extracted from pre-computed batch arrays
    lam_i_local       = lam_self(idx_i_mat(:, row));
    lam_p2_local      = lam_self(idx_j_mat(:, row));
    pair_qv_local     = qv_mat(:, row);                         % (2 x 1)
    tau_peanut_pair   = [tau_proj_top(:, row), tau_proj_bot(:, row)]; % N_c x 2

    use_dense_u_pair = use_ucorr_colloc ...
        && isfield(meta,'Ucross_colloc_actual') && ~isempty(meta.Ucross_colloc_actual) ...
        && isfield(meta,'Ec_colloc_actual')     && ~isempty(meta.Ec_colloc_actual);

    if use_dense_u_pair
        rhs_actual = [lam_i_local; lam_p2_local];
        u_fine     = meta.Ucross_colloc_actual * rhs_actual;
        if project_charge
            u_fine = u_fine - meta.Lr_colloc_actual * pair_qv_local;
        end
        u_peanut = meta.Ec_colloc_actual * tau_peanut_pair(:);
    else
        % Direct summation: recreate the close-pair interaction rhs
        u_fine_i  = -lapSLPfield(q(p2)+rbase_in_c, rcheck_out(block_i),  lam_p2_local, false);
        u_fine_p2 = -lapSLPfield(q(i) +rbase_in_c, rcheck_out(block_p2), lam_i_local,  false);
        if project_charge
            u_fine_i  = u_fine_i  - pair_qv_local(1);
            u_fine_p2 = u_fine_p2 - pair_qv_local(2);
        end
        u_fine   = [u_fine_i; u_fine_p2];
        u_peanut = lapSLPfield(rin_pair_c, rout_pair, tau_peanut_pair(:), false);
    end

    u_pair = u_fine - u_peanut;
    pair_idx = [block_i block_p2]';
    u_corr(pair_idx) = u_corr(pair_idx) + u_pair;
end

end % transform_peanut_laplace_fast


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function [top_out, bot_out] = batchRotatePair(top_in, bot_in, metas, ...
        rot_rows, N_c, is_real, phase_field)
%BATCHROTATEPAIR  Apply cached phase rotations to the two (N_c x n_pairs)
%  halves of a pair-input/output matrix.  Only non-identity pairs
%  (rot_rows) are touched; identity pairs pass through unchanged.
%
%  Supported phase types (auto-detected from the first non-identity pair):
%    numeric vector  - plain FFT-mode phase of length N_c
%    struct 'fft'    - spec.phase is length N_c
%    struct 'shift'  - converted to an equivalent FFT phase via spec.rot
%    struct 'oversampled_fft' - spec.phase is length M = n_oversampled;
%                               batch zero-padded FFT rotation

top_out = top_in;
bot_out = bot_in;
if isempty(rot_rows)
    return
end

n_rot = numel(rot_rows);

% Inspect the phase-spec type from the first non-identity pair
ph1 = metas(rot_rows(1)).(phase_field);

if isnumeric(ph1)
    % Plain numeric phase vector (length N_c) -- produced when
    % getUniformCircleRotationPhase is called directly.
    phases = zeros(N_c, n_rot);
    for k = 1:n_rot
        phases(:,k) = metas(rot_rows(k)).(phase_field)(:);
    end
    top_out(:, rot_rows) = applyPhaseRotation(top_in(:, rot_rows), phases, is_real);
    bot_out(:, rot_rows) = applyPhaseRotation(bot_in(:, rot_rows), phases, is_real);

elseif isstruct(ph1)
    switch ph1.mode

        case {'fft', 'shift'}
            % 'fft':   spec.phase is the N_c-length FFT phase vector.
            % 'shift': an exact circshift alias; convert to equivalent
            %          FFT phase for batched application (avoids grouping
            %          by shift_steps and handles all shifts uniformly).
            phases = zeros(N_c, n_rot);
            for k = 1:n_rot
                ph = metas(rot_rows(k)).(phase_field);
                if strcmp(ph.mode, 'fft')
                    phases(:,k) = ph.phase(:);
                else  % 'shift': derive the FFT phase from the rotation
                    phases(:,k) = getUniformCircleRotationPhase(N_c, ph.rot);
                end
            end
            top_out(:, rot_rows) = applyPhaseRotation(top_in(:, rot_rows), phases, is_real);
            bot_out(:, rot_rows) = applyPhaseRotation(bot_in(:, rot_rows), phases, is_real);

        case 'oversampled_fft'
            % Batch zero-padded FFT: zero-pad N_c spectrum to M, multiply
            % by M-length phase, IFFT, downsample back to N_c.  All n_rot
            % pairs are processed in a single pair of fft/ifft calls.
            M      = ph1.n_oversampled;
            phases = zeros(M, n_rot);
            for k = 1:n_rot
                phases(:,k) = metas(rot_rows(k)).(phase_field).phase(:);
            end
            top_out(:, rot_rows) = applyOversampledPhaseRotation( ...
                top_in(:, rot_rows), phases, N_c, M, is_real);
            bot_out(:, rot_rows) = applyOversampledPhaseRotation( ...
                bot_in(:, rot_rows), phases, N_c, M, is_real);

        otherwise
            % Unknown mode: serial fallback via rotateUniformCircleData
            for k = 1:n_rot
                row = rot_rows(k);
                top_out(:, row) = rotateUniformCircleData( ...
                    top_in(:, row), [], metas(row).(phase_field));
                bot_out(:, row) = rotateUniformCircleData( ...
                    bot_in(:, row), [], metas(row).(phase_field));
            end
    end
end
end % batchRotatePair

% -------------------------------------------------------------------------

function out = applyPhaseRotation(block, phases, force_real)
%APPLYPHASEROTATION  FFT-based phase rotation on a (N_c x n_rot) block.
%  block  : (N_c x n_rot) data to rotate
%  phases : (N_c x n_rot) precomputed Fourier phase vectors
%  out(:,k) = ifft(fft(block(:,k)) .* phases(:,k))
out = ifft(fft(block, [], 1) .* phases, [], 1);
if force_real
    out = real(out);
end
end % applyPhaseRotation

% -------------------------------------------------------------------------

function out = applyOversampledPhaseRotation(block, phases, N, M, force_real)
%APPLYOVERSAMPLEDPHASEROTATION  Batch zero-padded FFT phase rotation.
%  Replicates the per-column logic of rotateUniformCircleData 'oversampled_fft'
%  across all columns simultaneously.
%
%  block  : (N x n_rot) input data (real or complex)
%  phases : (M x n_rot) precomputed Fourier phase vectors at M = oversample*N
%  N, M   : fine and oversampled grid sizes
%  out    : (N x n_rot) rotated output

n_rot = size(block, 2);
half  = floor(N / 2);

% FFT of all columns at once
F   = fft(block, [], 1);     % (N x n_rot)

% Zero-pad spectrum from N to M  (replicating rotateUniformCircleData logic)
Fup = zeros(M, n_rot, 'like', F);
Fup(1:half+1, :)         = F(1:half+1, :);
Fup(M-N+half+2:M, :)     = F(half+2:N, :);
Fup = Fup * (M / N);

% Multiply by oversampled phase and IFFT
vals = ifft(Fup .* phases, [], 1);  % (M x n_rot)

% Downsample back to N points
stride = M / N;
out    = vals(1:stride:end, :);     % (N x n_rot)

if force_real
    out = real(out);
end
end % applyOversampledPhaseRotation
