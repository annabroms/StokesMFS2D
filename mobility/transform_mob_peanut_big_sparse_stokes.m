function [lam_c_x, lam_c_y, lambda_self_blocks, u_corr,...
    lam_c_nonpx, lam_c_nonpy] =...
    transform_mob_peanut_big_sparse_stokes(tau, geom, basis)
%TRANSFORM_MOB_PEANUT_BIG_SPARSE_STOKES Apply prebuilt sparse pair maps.
%
% This is the fast solve-grid, cmap/self-correct mobility peanut transform
% used by matvec_peanut_big_sparse. The first four outputs match the
% resistance big-sparse transform: coarse x/y sources, projected one-body
% sources in block form, and the close-pair velocity correction.
% Boundary postprocessing can additionally request the nonprojected coarse
% one-body sources.
%
% Implementation notes:
%   * lambda_self_blocks slices extracted once; all reshape outputs derived
%      from those two slices without re-indexing the block matrix.
%   *  applyDensePairProjector uses pagemtimes (R2020b+).
%   *  scatterSourceCorrections uses a single accumarray with pre-stacked
%      values.

if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('transform_mob_peanut_big_sparse_stokes:MissingBigSparse',...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck, geom.rvec_out)
    error('transform_mob_peanut_big_sparse_stokes:UnsupportedGrid',...
        'The big sparse transform is only supported on the solve grid.');
end

opt = geom.opt;
U   = basis.U;
Y   = basis.Y;
Lc  = basis.Lc;

P        = opt.P;
N_c      = opt.N_c;
PM       = round(P*N_c*opt.a_c);
N_large  = PM / P;
n_coarse = P*N_c;

%% 1-body: recover sources from boundary unknowns
% Apply U{1} then Y{1} to all particle columns at once (batched DGEMM).
tau_x      = reshape(tau(1:PM),      N_large, P);
tau_y      = reshape(tau(PM+1:2*PM), N_large, P);
tau_blocks = [tau_x; tau_y];

tau_mapped_blocks = Y{1} * (U{1} * tau_blocks);   % (2*N_c) x P

% Project: tau_proj = (I - Lc) * tau_mapped
tau_proj_blocks  = tau_mapped_blocks - Lc * tau_mapped_blocks;
lambda_self_blocks = tau_proj_blocks;

%% Extract x/y slices once; derive all four outputs from them
% index into tau_mapped_blocks and tau_proj_blocks once
% each, then reuse the column vectors for lambda_self assembly.
mapped_x = reshape(tau_mapped_blocks(1:N_c,:),       [], 1);  % (n_coarse x 1)
mapped_y = reshape(tau_mapped_blocks(N_c+1:2*N_c,:), [], 1);

proj_x   = reshape(tau_proj_blocks(1:N_c,:),         [], 1);
proj_y   = reshape(tau_proj_blocks(N_c+1:2*N_c,:),   [], 1);

lam_c_nonpx = mapped_x;
lam_c_nonpy = mapped_y;
lam_c_x     = proj_x;
lam_c_y     = proj_y;

% Assemble lambda_self from already-extracted vectors (no second indexing
% into lambda_self_blocks).
lambda_self = [proj_x; proj_y];   % (2*n_coarse x 1)

%% Source corrections
pair_proj = [];
use_sparse_map_coarse = getConfiguredSparseMapCoarse(opt, basis.big_sparse);

if use_sparse_map_coarse
    if isempty(basis.big_sparse.M_source_corr)
        error('transform_mob_peanut_big_sparse_stokes:MissingSourceCorr',...
            'opt.mob_sparse_map_coarse=1 requires M_source_corr.');
    end
    corr = basis.big_sparse.M_source_corr * lambda_self;
else
    if isempty(basis.big_sparse.M_pair_nonp)
        error('transform_mob_peanut_big_sparse_stokes:MissingPairMap',...
            'opt.mob_sparse_map_coarse=0 requires M_pair_nonp.');
    end
    if ~isfield(basis.big_sparse,'P_pair') || isempty(basis.big_sparse.P_pair)
        error('transform_mob_peanut_big_sparse_stokes:MissingProjector',...
            'opt.mob_sparse_map_coarse=0 requires basis.big_sparse.P_pair.');
    end
    if ~isfield(basis.big_sparse,'source_scatter_rows') ||...
            isempty(basis.big_sparse.source_scatter_rows)
        error('transform_mob_peanut_big_sparse_stokes:MissingSourceRows',...
            'opt.mob_sparse_map_coarse=0 requires source_scatter_rows.');
    end

    % Legacy source-correction algebra:
    %   pair_nonp = C_nonp * lambda_self
    %   pair_proj = P_pair * pair_nonp   (per-pair block projection)
    %   corr      = scatter([pair_proj; pair_nonp])
    pair_nonp = basis.big_sparse.M_pair_nonp * lambda_self;
    pair_proj = applyDensePairProjector(pair_nonp, basis.big_sparse.P_pair);
    corr      = scatterSourceCorrections(...
        basis.big_sparse.source_scatter_rows,...
        pair_proj, pair_nonp, 4*N_c, n_coarse);
end

lam_c_x = lam_c_x + corr(1:n_coarse);
lam_c_y = lam_c_y + corr(n_coarse+1:2*n_coarse);

%% Velocity correction
use_direct_u_corr = getConfiguredDirectUCorr(opt, basis.big_sparse);

if use_direct_u_corr
    if ~isfield(basis.big_sparse,'M_u_corr') ||...
            isempty(basis.big_sparse.M_u_corr)
        error('transform_mob_peanut_big_sparse_stokes:MissingUCorr',...
            'opt.big_sparse_direct_u_corr=1 requires M_u_corr.');
    end
    u_corr = basis.big_sparse.M_u_corr * lambda_self;
else
    if ~isfield(basis.big_sparse,'M_u_cross') ||...
            ~isfield(basis.big_sparse,'M_u_peanut') ||...
            isempty(basis.big_sparse.M_u_cross) ||...
            isempty(basis.big_sparse.M_u_peanut)
        error('transform_mob_peanut_big_sparse_stokes:MissingFactoredUCorr',...
            ['opt.big_sparse_direct_u_corr=0 requires M_u_cross ',...
             'and M_u_peanut.']);
    end
    % pair_proj is already computed in the sparse_map_coarse=false branch;
    % only compute it here when sparse_map_coarse=true.
    if isempty(pair_proj)
        if ~isfield(basis.big_sparse,'M_pair_proj') ||...
                isempty(basis.big_sparse.M_pair_proj)
            error('transform_mob_peanut_big_sparse_stokes:MissingPairProj',...
                ['opt.mob_sparse_map_coarse=1 and opt.big_sparse_direct_u_corr=0 ',...
                 'require M_pair_proj.']);
        end
        pair_proj = basis.big_sparse.M_pair_proj * lambda_self;
    end
    u_corr = basis.big_sparse.M_u_cross  * lambda_self -...
             basis.big_sparse.M_u_peanut * pair_proj;
end

end % transform_mob_peanut_big_sparse_stokes

% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function tf = getConfiguredSparseMapCoarse(opt, big_sparse)
tf = logical(getOptField(opt,'mob_sparse_map_coarse',...
    getOptField(opt,'sparse_map_coarse',...
    getBigSparsePlanField(big_sparse,'sparse_map_coarse',true))));
end

function tf = getConfiguredDirectUCorr(opt, big_sparse)
tf = logical(getOptField(opt,'big_sparse_direct_u_corr',...
    getBigSparsePlanField(big_sparse,'direct_u_corr',true)));
end

function value = getBigSparsePlanField(big_sparse, name, default_value)
if isfield(big_sparse,'matrix_plan') && isfield(big_sparse.matrix_plan, name)
    value = big_sparse.matrix_plan.(name);
else
    value = default_value;
end
end

% -------------------------------------------------------------------------

function pair_proj = applyDensePairProjector(pair_nonp, P_pair)
%APPLYDENSEPAIRPROJECTOR  Per-pair block projection via pagemtimes.
%
%  pair_nonp : (pair_rows*n_pairs x 1)  – flat column vector
%  P_pair    : (pair_rows x pair_rows)  – dense projector applied per pair
%  pair_proj : (pair_rows*n_pairs x 1)  – flat column vector

if isempty(pair_nonp)
    pair_proj = pair_nonp;
    return
end

pair_rows = size(P_pair, 1);
n_pairs   = numel(pair_nonp) / pair_rows;

% Reshape to (pair_rows x n_pairs), multiply, flatten.
% For n_pairs == 1 this is a plain matrix-vector product.
% For n_pairs > 1, pagemtimes broadcasts P_pair across all pair columns.
% Expand P_pair to (pair_rows x pair_rows x 1) so pagemtimes broadcasts.
pair_proj = pagemtimes(P_pair, reshape(pair_nonp, pair_rows, n_pairs));
pair_proj = pair_proj(:);
end

% -------------------------------------------------------------------------

function corr = scatterSourceCorrections(rows, pair_proj, pair_nonp,...
    pair_rows, n_coarse)
%SCATTERSOURCECORRECTIONS Scatter [pair_proj; pair_nonp] into a correction
% vector of length 4*n_coarse using accumarray.
%
% Stack pair_proj and pair_nonp into a single flat value
% vector in one operation, then call accumarray once.

n_source_rows = 4 * n_coarse;

if isempty(rows)
    corr = zeros(n_source_rows, 1);
    return
end

n_pairs = numel(pair_nonp) / pair_rows;

% Stack projected and nonprojected values in pair-interleaved order that
% matches the row index layout in `rows`, then scatter in one accumarray.
vals = [reshape(pair_proj, pair_rows, n_pairs);...
        reshape(pair_nonp, pair_rows, n_pairs)];   % (2*pair_rows x n_pairs)
corr = accumarray(rows(:), vals(:), [n_source_rows, 1]);
end
