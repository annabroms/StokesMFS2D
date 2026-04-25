function [lam_c, lam_self, u_corr, pair_qv_nonp,...
    lam_c_nonp, lam_self_nonp] =...
    transform_lap_peanut_big_sparse(tau, geom, basis)
%TRANSFORM_LAP_PEANUT_BIG_SPARSE Apply prebuilt sparse pair maps.
%
% Outputs:
%   lam_c         - projected coarse sources (self + pair correction)
%   lam_self      - projected one-body coarse sources only
%   u_corr        - close-pair velocity correction on the solve grid
%   pair_qv_nonp  - scattered non-projected pair charge/voltage
%   lam_c_nonp    - non-projected coarse sources (self + pair correction)
%   lam_self_nonp - non-projected one-body coarse sources only

if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('transform_lap_peanut_big_sparse:MissingBigSparse',...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck, geom.rvec_out)
    error('transform_lap_peanut_big_sparse:UnsupportedGrid',...
        'The Laplace big-sparse transform is only supported on the solve grid.');
end

opt      = geom.opt;
rvec_out = geom.rvec_out;
P        = opt.P;
N_c      = opt.N_c;
N_large  = numel(rvec_out) / P;
n_coarse = P * N_c;

project_charge = logical(getOptField(opt,'project_charge',false));

%% 1-body basis
% All one-body blocks are identical, so batch them as particle columns
% while keeping the U{1} then Y{1} ordering for a backward-stable apply.
tau_blocks    = reshape(tau(1:P*N_large), N_large, P);
lam_blocks_nonp = basis.Y{1} * (basis.U{1} * tau_blocks);  % (N_c x P)

% Charge projection: subtract per-column mean if requested
if project_charge
    lam_blocks = lam_blocks_nonp - mean(lam_blocks_nonp, 1);
else
    lam_blocks = lam_blocks_nonp;
end

% lam_self / lam_self_nonp capture the one-body result BEFORE pair
% corrections are added. lam_c / lam_c_nonp will be updated below.
lam_self      = lam_blocks(:);          % (n_coarse x 1)
lam_self_nonp = lam_blocks_nonp(:);

%% 2-body basis
% Apply pair maps to the projected one-body sources, project the pair
% sources, then scatter corrections back to coarse sources and velocity.
[pair_nonp, pair_qv] = applyPairMaps(basis.big_sparse, lam_self);
pair_proj = applyDensePairProjector(pair_nonp, basis.big_sparse.P_pair);

corr = scatterSourceCorrections(...
    basis.big_sparse.source_scatter_rows,...
    pair_proj, pair_nonp, 2*N_c, n_coarse);

lam_c      = lam_self      + corr(1:n_coarse);
lam_c_nonp = lam_self_nonp + corr(n_coarse+1:2*n_coarse);

pair_qv_nonp = scatterPairQv(basis.big_sparse.qv_scatter_rows, pair_qv, P);

u_corr = basis.big_sparse.M_u_cross  * lam_self -...
         basis.big_sparse.M_u_peanut * pair_proj;

% Optional charge/voltage correction (field present and non-empty)
if ~isempty(basis.big_sparse.M_u_qv)
    u_corr = u_corr - basis.big_sparse.M_u_qv * pair_qv;
end

end % transform_lap_peanut_big_sparse

% =========================================================================
%  PAIR MAP APPLICATION
% =========================================================================

function [pair_nonp, pair_qv] = applyPairMaps(big_sparse, lam_self)
%APPLYPAIRMAPS Apply the pair source maps, with optional rotation reuse.
if isfield(big_sparse,'pair_map_mode') &&...
        strcmp(big_sparse.pair_map_mode,'rotated_factored')
    pair_rhs      = rotatePairInputs(lam_self, big_sparse);
    pair_nonp     = rotatePairOutputs(...
        big_sparse.M_pair_nonp_canon * pair_rhs, big_sparse);
    pair_qv       = big_sparse.M_pair_qv_canon * pair_rhs;
else
    pair_nonp = big_sparse.M_pair_nonp * lam_self;
    pair_qv   = big_sparse.M_pair_qv  * lam_self;
end
end

% -------------------------------------------------------------------------

function pair_rhs = rotatePairInputs(lam_self, big_sparse)
%ROTATEPAIRINPUTS Gather and batch-rotate one-body sources for all pairs.
%
% Uses precomputed phase_c_mat (N_c x n_pairs) from big_sparse when
% available, avoiding all per-matvec phase extraction cost.
% Falls back to a stacking loop if the cache is absent, and to serial
% rotateUniformCircleData calls if phase_c is a struct (spec mode).

N_c     = big_sparse.N_c;
pairs   = big_sparse.pairs;
meta    = big_sparse.pair_rotation_meta;

% --- Vectorised gather: (2*N_c) x n_pairs ---
local = (1:N_c)';
idx_i = bsxfun(@plus, (pairs(:,1)-1)'*N_c, local);   % N_c x n_pairs
idx_j = bsxfun(@plus, (pairs(:,2)-1)'*N_c, local);
pair_rhs_mat = [lam_self(idx_i); lam_self(idx_j)];   % (2*N_c) x n_pairs

% --- Vectorised identity check ---
rots      = [meta.rot];
needs_rot = abs(rots - 1) > 100 * eps(max(1, abs(rots)));
rot_rows  = find(needs_rot);

if ~isempty(rot_rows)
    % Resolve phase matrix: cached (O(1) slice) > stacking loop > serial
    has_cache = isfield(big_sparse,'phase_c_mat') &&...
        ~isempty(big_sparse.phase_c_mat);

    if has_cache
        phases    = big_sparse.phase_c_mat(:, rot_rows);   % N_c x n_rot
        use_batch = true;
    else
        use_batch = isnumeric(meta(rot_rows(1)).phase_c);
        if use_batch
            n_rot  = numel(rot_rows);
            phases = zeros(N_c, n_rot);
            for k = 1:n_rot
                phases(:,k) = meta(rot_rows(k)).phase_c(:);
            end
        end
    end

    if use_batch
        is_real = isreal(lam_self);
        pair_rhs_mat(1:N_c,     rot_rows) = applyPhaseRotation(...
            pair_rhs_mat(1:N_c,     rot_rows), phases, is_real);
        pair_rhs_mat(N_c+1:end, rot_rows) = applyPhaseRotation(...
            pair_rhs_mat(N_c+1:end, rot_rows), phases, is_real);
    else
        % Struct-spec fallback: serial call per non-identity pair
        for row = rot_rows
            rhs     = reshape(pair_rhs_mat(:,row), N_c, 2);
            rotated = rotateUniformCircleData(rhs, [], meta(row).phase_c);
            pair_rhs_mat(:,row) = rotated(:);
        end
    end
end

pair_rhs = pair_rhs_mat(:);
end

% -------------------------------------------------------------------------

function pair_nonp = rotatePairOutputs(pair_nonp_canon, big_sparse)
%ROTATEPAIROUTPUTS Batch inverse-rotate pair outputs back to lab frame.
%
% Uses precomputed phase_c_inv_mat (N_c x n_pairs) from big_sparse when
% available. Falls back to stacking loop or serial calls.

N_c     = big_sparse.N_c;
meta    = big_sparse.pair_rotation_meta;
n_pairs = numel(pair_nonp_canon) / (2*N_c);

% --- Vectorised identity check ---
rots      = [meta.rot];
needs_rot = abs(rots - 1) > 100 * eps(max(1, abs(rots)));
rot_rows  = find(needs_rot);

% Identity pairs are already correct: initialise output as a copy and
% only overwrite non-identity pairs.
pair_nonp = pair_nonp_canon;

if ~isempty(rot_rows)
    has_cache = isfield(big_sparse,'phase_c_inv_mat') &&...
        ~isempty(big_sparse.phase_c_inv_mat);

    if has_cache
        phases    = big_sparse.phase_c_inv_mat(:, rot_rows);  % N_c x n_rot
        use_batch = true;
    else
        use_batch = isnumeric(meta(rot_rows(1)).phase_c_inv);
        if use_batch
            n_rot  = numel(rot_rows);
            phases = zeros(N_c, n_rot);
            for k = 1:n_rot
                phases(:,k) = meta(rot_rows(k)).phase_c_inv(:);
            end
        end
    end

    canon_mat = reshape(pair_nonp_canon, 2*N_c, n_pairs);

    if use_batch
        is_real = isreal(pair_nonp_canon);
        out_mat = reshape(pair_nonp, 2*N_c, n_pairs);
        out_mat(1:N_c,     rot_rows) = applyPhaseRotation(...
            canon_mat(1:N_c,     rot_rows), phases, is_real);
        out_mat(N_c+1:end, rot_rows) = applyPhaseRotation(...
            canon_mat(N_c+1:end, rot_rows), phases, is_real);
        pair_nonp = out_mat(:);
    else
        % Struct-spec fallback: serial call per non-identity pair
        for row = rot_rows
            pair_idx            = (row-1)*2*N_c+1 : row*2*N_c;
            tau_pair            = reshape(canon_mat(:,row), N_c, 2);
            rotated             = rotateUniformCircleData(...
                tau_pair, [], meta(row).phase_c_inv);
            pair_nonp(pair_idx) = rotated(:);
        end
    end
end
end

% -------------------------------------------------------------------------

function tf = isIdentityRotation(rot)
tf = abs(rot - 1) <= 100 * eps(max(1, abs(rot)));
end

% =========================================================================
%  PAIR PROJECTOR
% =========================================================================

function pair_proj = applyDensePairProjector(pair_nonp, P_pair)
%APPLYDENSEPAIRPROJECTOR Per-pair block projection via pagemtimes (R2020b+).
%
%  pair_nonp : (pair_rows*n_pairs x 1)
%  P_pair    : (pair_rows x pair_rows)
%  pair_proj : (pair_rows*n_pairs x 1)

if isempty(pair_nonp)
    pair_proj = pair_nonp;
    return
end

pair_rows = size(P_pair, 1);
n_pairs   = numel(pair_nonp) / pair_rows;

% pagemtimes broadcasts P_pair across all pair columns in one call.
pair_proj = pagemtimes(P_pair, reshape(pair_nonp, pair_rows, n_pairs));
pair_proj = pair_proj(:);
end

% =========================================================================
%  SCATTER HELPERS
% =========================================================================

function corr = scatterSourceCorrections(rows, pair_proj, pair_nonp,...
    pair_rows, n_coarse)
%SCATTERSOURCECORRECTIONS Scatter [pair_proj; pair_nonp] into a correction
% vector of length 2*n_coarse using a single accumarray call.

if isempty(rows)
    corr = zeros(2*n_coarse, 1);
    return
end

n_pairs = numel(pair_nonp) / pair_rows;

% Stack projected and non-projected values in the order that matches rows,
% then scatter in one accumarray call.
vals = [reshape(pair_proj, pair_rows, n_pairs);
        reshape(pair_nonp, pair_rows, n_pairs)];   % (2*pair_rows x n_pairs)
corr = accumarray(rows(:), vals(:), [2*n_coarse, 1]);
end

% -------------------------------------------------------------------------

function pair_qv_nonp = scatterPairQv(rows, pair_qv, P)
%SCATTERPAIRQV Scatter pair charge/voltage values to per-particle totals.
if isempty(rows)
    pair_qv_nonp = zeros(P, 1);
else
    pair_qv_nonp = accumarray(rows(:), pair_qv(:), [P, 1]);
end
end

function out = applyPhaseRotation(block, phases, force_real)
%APPLYPHASEROTATION FFT-based phase rotation on a (N_c x n_rot) block.
%
%  block      : (N_c x n_rot) data to rotate
%  phases     : (N_c x n_rot) precomputed Fourier phase vectors
%  force_real : scalar logical — strip imaginary part if input was real
%
%  out(:,k) = ifft(fft(block(:,k)).* phases(:,k))
%
%  Single fft/ifft call regardless of how many columns are present.

out = ifft(fft(block, [], 1).* phases, [], 1);
if force_real
    out = real(out);
end
end





% function [lam_c,lam_self,u_corr,pair_qv_nonp, ...
%     lam_c_nonp,lam_self_nonp] = ...
%     transform_lap_peanut_big_sparse(tau,geom,basis)
% %TRANSFORM_LAP_PEANUT_BIG_SPARSE Apply prebuilt sparse pair maps.
% 
% if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
%     error('transform_lap_peanut_big_sparse:MissingBigSparse', ...
%         'basis.big_sparse is required for opt.use_big_sparse=1.');
% end
% if ~isequal(geom.rcheck,geom.rvec_out)
%     error('transform_lap_peanut_big_sparse:UnsupportedGrid', ...
%         'The Laplace big-sparse transform is only supported on the solve grid.');
% end
% 
% opt = geom.opt;
% rvec_out = geom.rvec_out;
% U = basis.U;
% Y = basis.Y;
% 
% P = opt.P;
% N_c = opt.N_c;
% N_large = numel(rvec_out)/P;
% n_coarse = P*N_c;
% project_charge = logical(getOptField(opt,'project_charge',false));
% 
% % 1-body basis. All one-body blocks are identical, so batch them as
% % particle columns while keeping the U{1} then Y{1} ordering for a backward stable apply.
% tau_blocks = reshape(tau(1:P*N_large),N_large,P);
% lam_blocks_nonp = Y{1}*(U{1}*tau_blocks);
% if project_charge
%     lam_blocks = lam_blocks_nonp - mean(lam_blocks_nonp,1);
% else
%     lam_blocks = lam_blocks_nonp;
% end
% lam_c = reshape(lam_blocks,[],1);
% lam_c_nonp = reshape(lam_blocks_nonp,[],1);
% 
% lam_self = lam_c;
% lam_self_nonp = lam_c_nonp;
% 
% % 2-body basis. Apply the pair maps to the projected one-body sources to get the close-pair contributions, then scatter them back to the coarse sources and velocity correction.
% [pair_nonp,pair_qv] = applyPairMaps(basis.big_sparse,lam_self);
% pair_proj = applyDensePairProjector(pair_nonp,basis.big_sparse.P_pair);
% corr = scatterSourceCorrections(basis.big_sparse.source_scatter_rows, ...
%     pair_proj,pair_nonp,2*N_c,n_coarse);
% 
% lam_c = lam_self + corr(1:n_coarse);
% lam_c_nonp = lam_self_nonp + corr(n_coarse+1:2*n_coarse);
% 
% pair_qv_nonp = scatterPairQv(basis.big_sparse.qv_scatter_rows,pair_qv,P);
% u_corr = basis.big_sparse.M_u_cross*lam_self - ...
%     basis.big_sparse.M_u_peanut*pair_proj;
% if isfield(basis.big_sparse,'M_u_qv') && ~isempty(basis.big_sparse.M_u_qv)
%     u_corr = u_corr - basis.big_sparse.M_u_qv*pair_qv;
% end
% 
% 
% end
% 
% function [pair_nonp,pair_qv] = applyPairMaps(big_sparse,lam_self)
% if isfield(big_sparse,'pair_map_mode') && ...
%         strcmp(big_sparse.pair_map_mode,'rotated_factored')
%     pair_rhs = rotatePairInputs(lam_self,big_sparse);
%     pair_nonp_canon = big_sparse.M_pair_nonp_canon*pair_rhs;
%     pair_nonp = rotatePairOutputs(pair_nonp_canon,big_sparse);
%     pair_qv = big_sparse.M_pair_qv_canon*pair_rhs;
% else
%     pair_nonp = big_sparse.M_pair_nonp*lam_self;
%     pair_qv = big_sparse.M_pair_qv*lam_self;
% end
% end
% 
% function pair_rhs = rotatePairInputs(lam_self,big_sparse)
% N_c = big_sparse.N_c;
% pairs = big_sparse.pairs;
% meta = big_sparse.pair_rotation_meta;
% n_pairs = size(pairs,1);
% pair_rhs = zeros(n_pairs*2*N_c,1);
% 
% for row = 1:n_pairs
%     i = pairs(row,1);
%     j = pairs(row,2);
%     idx_i = (i-1)*N_c+1:i*N_c;
%     idx_j = (j-1)*N_c+1:j*N_c;
%     rhs = [lam_self(idx_i) lam_self(idx_j)];
%     if ~isIdentityRotation(meta(row).rot)
%         rhs = rotateUniformCircleData(rhs,[],meta(row).phase_c);
%     end
%     pair_idx = (row-1)*2*N_c+1:row*2*N_c;
%     pair_rhs(pair_idx) = rhs(:);
% end
% end
% 
% function pair_nonp = rotatePairOutputs(pair_nonp_canon,big_sparse)
% N_c = big_sparse.N_c;
% meta = big_sparse.pair_rotation_meta;
% n_pairs = numel(pair_nonp_canon)/(2*N_c);
% pair_nonp = zeros(size(pair_nonp_canon));
% 
% for row = 1:n_pairs
%     pair_idx = (row-1)*2*N_c+1:row*2*N_c;
%     tau_pair = reshape(pair_nonp_canon(pair_idx),N_c,2);
%     if ~isIdentityRotation(meta(row).rot)
%         tau_pair = rotateUniformCircleData(tau_pair,[], ...
%             meta(row).phase_c_inv);
%     end
%     pair_nonp(pair_idx) = tau_pair(:);
% end
% end
% 
% function tf = isIdentityRotation(rot)
% tf = abs(rot - 1) <= 100*eps(max(1,abs(rot)));
% end
% 
% function pair_proj = applyDensePairProjector(pair_nonp,P_pair)
% if isempty(pair_nonp)
%     pair_proj = pair_nonp;
%     return
% end
% pair_rows = size(P_pair,1);
% n_pairs = numel(pair_nonp)/pair_rows;
% pair_proj = P_pair*reshape(pair_nonp,pair_rows,n_pairs);
% pair_proj = pair_proj(:);
% end
% 
% function corr = scatterSourceCorrections(rows,pair_proj,pair_nonp, ...
%     pair_rows,n_coarse)
% if isempty(rows)
%     corr = zeros(2*n_coarse,1);
%     return
% end
% n_pairs = numel(pair_nonp)/pair_rows;
% pair_proj = reshape(pair_proj,pair_rows,n_pairs);
% pair_nonp = reshape(pair_nonp,pair_rows,n_pairs);
% vals = reshape([pair_proj; pair_nonp],[],1);
% corr = accumarray(rows,vals,[2*n_coarse,1]);
% end
% 
% function pair_qv_nonp = scatterPairQv(rows,pair_qv,P)
% if isempty(rows)
%     pair_qv_nonp = zeros(P,1);
% else
%     pair_qv_nonp = accumarray(rows,pair_qv(:),[P,1]);
% end
% end
% 
% function lam_out = projectChargeMode(lam_in,project_charge)
% if ~project_charge || isempty(lam_in)
%     lam_out = lam_in;
%     return
% end
% n = numel(lam_in);
% lam_out = lam_in - (sum(lam_in)/n);
% end
