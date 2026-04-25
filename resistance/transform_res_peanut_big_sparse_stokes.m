function [lam_c_x, lam_c_y, lambda_self_blocks, u_corr] =...
    transform_res_peanut_big_sparse_stokes(tau,geom,basis)
%TRANSFORM_RES_PEANUT_BIG_SPARSE_STOKES Apply one-body pseudoinverse plus
% solve-grid sparse pair maps.
%
% Optimizations applied (Apr 2026):
%   1. Nii applied once outside the chunk loop (2 DGEMM calls total).
%   2. applyDensePageBlocks uses pagemtimes (R2020b+) instead of bsxfun.
%   3. lambda_body column-gathers precomputed outside the chunk loop.
%   4. addSourceCorrectionsInPairOrder / addVectorCorrectionsInPairOrder
%      replaced with vectorised accumarray scatter-adds (no loops).

if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('transform_res_peanut_big_sparse_stokes:MissingBigSparse',...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('transform_res_peanut_big_sparse_stokes:UnsupportedGrid',...
        'The resistance big sparse transform is only supported on the solve grid.');
end

opt      = geom.opt;
rvec_out = geom.rvec_out;

P        = opt.P;
N_c      = opt.N_c;
PM       = numel(rvec_out);
N_large  = PM/P;
n_coarse = P*N_c;

%% One-body pseudoinverse (batched over all P bodies)
tau_x        = reshape(tau(1:PM),       N_large, P);
tau_y        = reshape(tau(PM+1:2*PM),  N_large, P);
tau_blocks   = [tau_x; tau_y];

step_blocks          = basis.U{1} * tau_blocks;
lambda_self_blocks   = basis.Y{1} * step_blocks;

lam_self_x   = reshape(lambda_self_blocks(1:N_c,:),         [], 1);
lam_self_y   = reshape(lambda_self_blocks(N_c+1:2*N_c,:),   [], 1);
lambda_self  = [lam_self_x; lam_self_y];

%% Pair corrections through the fixed sparse map
pair_tau             = basis.big_sparse.M_pair * lambda_self;
use_sparse_map_coarse = getConfiguredSparseMapCoarse(opt, basis.big_sparse);

if use_sparse_map_coarse
    if ~isfield(basis.big_sparse,'M_source_corr') ||...
            isempty(basis.big_sparse.M_source_corr)
        error('transform_res_peanut_big_sparse_stokes:MissingSourceCorr',...
            'opt.res_sparse_map_coarse=1 requires M_source_corr.');
    end
    lambda_c = lambda_self + basis.big_sparse.M_source_corr * lambda_self;
else
    if ~isfield(basis.big_sparse,'source_scatter_rows') ||...
            isempty(basis.big_sparse.source_scatter_rows)
        error('transform_res_peanut_big_sparse_stokes:MissingSourceRows',...
            'opt.res_sparse_map_coarse=0 requires source_scatter_rows.');
    end
    % Vectorised scatter-add (no loop over pairs)
    lambda_c = lambda_self;
    rows     = basis.big_sparse.source_scatter_rows;
    lambda_c = lambda_c + accumarray(rows(:), pair_tau(:), size(lambda_c));
end

lam_c_x = lambda_c(1:n_coarse);
lam_c_y = lambda_c(n_coarse+1:2*n_coarse);

%% Velocity field correction
u_corr_mode    = getVelocityCorrectionMode(opt, basis.big_sparse);
use_small_u_corr = getConfiguredResSmallerUCorr(opt, basis.big_sparse);

switch u_corr_mode
    case 'combined'
        if isempty(pair_tau)
            pair_u_corr = zeros(0,1);
        elseif use_small_u_corr
            if ~isfield(basis.big_sparse,'u_corr_small') ||...
                    isempty(basis.big_sparse.u_corr_small)
                error('transform_res_peanut_big_sparse_stokes:MissingSmallUCorr',...
                    'opt.res_smaller_ucorr=1 requires u_corr_small.');
            end
            pair_u_corr = applySmallCombinedUCorr(...
                lambda_self, pair_tau, geom, basis);
        else
            if ~isfield(basis.big_sparse,'M_u_corr') ||...
                    isempty(basis.big_sparse.M_u_corr)
                error('transform_res_peanut_big_sparse_stokes:MissingUCorr',...
                    'opt.res_smaller_ucorr=0 requires M_u_corr.');
            end
            pair_u_corr = basis.big_sparse.M_u_corr * [lambda_self; pair_tau];
        end
    otherwise
        if ~isfield(basis.big_sparse,'M_u_cross') ||...
                ~isfield(basis.big_sparse,'M_u_peanut') ||...
                isempty(basis.big_sparse.M_u_cross) ||...
                isempty(basis.big_sparse.M_u_peanut)
            error('transform_res_peanut_big_sparse_stokes:MissingFactoredUCorr',...
                ['opt.res_big_sparse_u_corr_mode=''factored'' requires ',...
                 'M_u_cross and M_u_peanut.']);
        end
        pair_u_corr = basis.big_sparse.M_u_cross  * lambda_self -...
                      basis.big_sparse.M_u_peanut * pair_tau;
end

% Vectorised scatter-add for velocity corrections (no loop over pairs)
out  = zeros(2*numel(geom.rcheck), 1);
rows = basis.big_sparse.velocity_scatter_rows;
u_corr = out + accumarray(rows(:), pair_u_corr(:), size(out));

end % transform_res_peanut_big_sparse_stokes

% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function pair_u_corr = applySmallCombinedUCorr(lambda_self, pair_tau, geom, basis)
%APPLYSMALLCOMBINEDCORR Optimised small combined u-correction.
%
%  Key optimisations vs. original:
%    * basis.Nii applied to the full tau_i / tau_j matrices in two DGEMM
%      calls before the chunk loop, eliminating per-chunk Nii multiplies.
%    * Column gathers lambda_body(:,pairs(:,1/2)) precomputed outside loop.
%    * applyDensePageBlocks uses pagemtimes instead of bsxfun+sum.

opt      = geom.opt;
P        = opt.P;
N_c      = opt.N_c;
N_check  = numel(geom.rcheck) / P;
n_coarse = P * N_c;

pairs   = basis.big_sparse.pairs;
payload = basis.big_sparse.u_corr_small;
n_pairs = size(pairs, 1);

%% Reshape sources
lambda_body = [reshape(lambda_self(1:n_coarse),           N_c, P);...
               reshape(lambda_self(n_coarse+1:2*n_coarse), N_c, P)];

pair_tau_blocks = reshape(pair_tau, 4*N_c, n_pairs);
[tau_i, tau_j]  = pairOrderToBodyBlocks(pair_tau_blocks, N_c);

%% Optimisation 1 – hoist Nii outside the loop (2 DGEMM calls total)
Nii_tau_i = basis.Nii * tau_i;   % (2*N_c) x n_pairs
Nii_tau_j = basis.Nii * tau_j;   % (2*N_c) x n_pairs

%% Optimisation 2 – precompute column gathers
lambda_body_i = lambda_body(:, pairs(:,1));  % (2*N_c) x n_pairs
lambda_body_j = lambda_body(:, pairs(:,2));  % (2*N_c) x n_pairs

%% Chunked loop (memory-limited; only page-block contractions remain)
pair_u_blocks = zeros(4*N_check, n_pairs);
chunk_pairs   = max(1, payload.chunk_pairs);

for first = 1:chunk_pairs:n_pairs
    last = min(first + chunk_pairs - 1, n_pairs);
    idx  = first:last;

    rhs_ij = lambda_body_j(:,idx) + tau_j(:,idx);
    rhs_ji = lambda_body_i(:,idx) + tau_i(:,idx);

    u_i = -Nii_tau_i(:,idx) - applyDensePageBlocks(payload.cross_ij(:,:,idx), rhs_ij);
    u_j = -Nii_tau_j(:,idx) - applyDensePageBlocks(payload.cross_ji(:,:,idx), rhs_ji);

    pair_u_blocks(:,idx) = bodyBlocksToPairOrder(u_i, u_j, N_check);
end

pair_u_corr = pair_u_blocks(:);
end

% -------------------------------------------------------------------------

function tf = getConfiguredSparseMapCoarse(opt, big_sparse)
tf = logical(getOptField(opt,'res_sparse_map_coarse',...
    getOptField(opt,'sparse_map_coarse',...
    getBigSparsePlanField(big_sparse,'sparse_map_coarse',false))));
end

function tf = getConfiguredResSmallerUCorr(opt, big_sparse)
tf = logical(getOptField(opt,'res_smaller_ucorr',...
    getOptField(opt,'smaller_ucorr',...
    getBigSparsePlanField(big_sparse,'smaller_u_corr',false))));
end

function value = getBigSparsePlanField(big_sparse, name, default_value)
if isfield(big_sparse,'matrix_plan') && isfield(big_sparse.matrix_plan, name)
    value = big_sparse.matrix_plan.(name);
else
    value = default_value;
end
end

function mode = getVelocityCorrectionMode(opt, big_sparse)
if isfield(opt,'res_big_sparse_u_corr_mode')
    mode = resolveResBigSparseUCorrMode(opt);
elseif isfield(big_sparse,'matrix_plan') &&...
        isfield(big_sparse.matrix_plan,'u_corr_mode')
    mode = big_sparse.matrix_plan.u_corr_mode;
else
    mode = resolveResBigSparseUCorrMode(opt);
end
end

% -------------------------------------------------------------------------
%  Pair-order / body-block layout helpers
% -------------------------------------------------------------------------

function [body_i, body_j] = pairOrderToBodyBlocks(pair_blocks, n_rows)
body_i = [pair_blocks(1:n_rows,:);           pair_blocks(2*n_rows+1:3*n_rows,:)];
body_j = [pair_blocks(n_rows+1:2*n_rows,:);  pair_blocks(3*n_rows+1:4*n_rows,:)];
end

function pair_blocks = bodyBlocksToPairOrder(body_i, body_j, n_rows)
pair_blocks = [body_i(1:n_rows,:);         body_j(1:n_rows,:);...
               body_i(n_rows+1:2*n_rows,:); body_j(n_rows+1:2*n_rows,:)];
end

% -------------------------------------------------------------------------
%  Dense page-block matrix-vector product
%  Optimisation 3 – pagemtimes (R2020b+) replaces bsxfun+sum
% -------------------------------------------------------------------------

function y = applyDensePageBlocks(Apages, x)
%APPLYDENSEPAGEBLOCKS  Batched matrix-vector product  y(:,k) = Apages(:,:,k)*x(:,k).
%
%  Apages : (m x n x n_pages)
%  x      : (n x n_pages)
%  y      : (m x n_pages)

n_pages = size(Apages, 3);
if n_pages == 0
    y = zeros(size(Apages,1), 0);
    return
end

% Reshape x to (n x 1 x n_pages) so pagemtimes treats each page as a
% matrix-vector product, then squeeze the singleton dimension.
y = pagemtimes(Apages, reshape(x, size(x,1), 1, n_pages));
y = reshape(y, size(Apages,1), n_pages);
end
