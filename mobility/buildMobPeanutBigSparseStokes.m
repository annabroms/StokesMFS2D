function [big_sparse,stats] = buildMobPeanutBigSparseStokes(geom,basis)
%BUILDMOBPEANUTBIGSPARSESTOKES Build solve-grid sparse pair-correction
%maps. These are both to correct coarse sources and to correct the velocity 
% field on the close pairs themselves due to coarse sources
%
% Calling with no inputs runs a visual self-test on a small random_discs_mc
% packing at phi=0.65 and shows the sparse matrix structures with spy and,
% if requested, the actual magnitudes of the compressed matrices.
%
% This builder is deliberately actual-pair only (no canonical mappings for
% shared rotations). solve_mob_peanut_enhanced disables
% opt.reuse_pair_basis_by_sep when opt.use_big_sparse is true, so every
% close pair contributes the dense maps for its actual geometry.
%
% opt.mob_sparse_map_coarse and opt.big_sparse_direct_u_corr select the
% matrix set built here:
%   mob_sparse_map_coarse = 0:
%       Coarse corrections are computed using pair-based row structure,
%       projecting pair non-projected sources and scattering the results to
%       the particles within the matvec.
%   mob_sparse_map_coarse = 1:
%       Coarse corrections are computed by directly mapping projected
%       one-body coarse sources to projected coarse-source increments
%       (main difference: this is per particle, NOT per pair).
%   big_sparse_direct_u_corr = 1:
%       Velocity corrections are computed by directly mapping projected
%       one-body coarse sources to velocity corrections.
%   big_sparse_direct_u_corr = 0:
%       Velocity corrections are factored into cross and pair components.
%       When mob_sparse_map_coarse = 1, pair projection data is precomputed
%       to avoid runtime dense projection.
%   Rigid-body motion corrections are always built for postprocessing.
%
% Note: The fastest and least memory-heavy setup is big_sparse_direct_u_corr = 1
% and mob_sparse_map_coarse = 1. Setting big_sparse_direc_u_corr = 0 might be
% slighly better behaved numerically, but the difference seems very tiny.
%
% All the options described here are compared in the script
% apr24_mobility_sparse_map_coarse_compare.m in the experiments/ folder.

if nargin == 0
    test_buildMobPeanutBigSparseStokes;
    big_sparse = [];
    stats      = [];
    return
end

timer = tic;

opt      = geom.opt;
q        = geom.q(:);
pairs    = geom.pairs;
rcheck   = geom.rcheck;
N_c      = opt.N_c;
P        = numel(q);
n_pairs  = size(pairs,1);
n_coarse = P * N_c;
N_check  = numel(rcheck) / P;

build_mode = resolveMobBigSparseBuildMode(opt);
validateBigSparseInputs(geom, basis, N_check, build_mode);

%% Dimension constants (derived once, shared throughout)
n_source_cols  = 2 * n_coarse;
n_u_rows       = 2 * P * N_check;
n_u_cols       = 2 * n_coarse;
pair_block_rows = 4 * N_c;
pair_total_rows = n_pairs * pair_block_rows;
n_rbm_rows     = 3 * P;

%% RAM plan
ram_estimate      = estimateMobPeanutBigSparseRamStokes(P,N_c,N_check,n_pairs,opt);
plan              = ram_estimate.matrix_plan;
direct_u_corr     = plan.direct_u_corr;
sparse_map_coarse = plan.sparse_map_coarse;

%% Pair projection matrix
if n_pairs > 0
    P_pair = resolvePairProjectionMatrix(geom, basis);
else
    P_pair = [];
end

%% Stats initialisation & build-memory guard
stats           = initBigSparseStats(n_pairs, N_c, N_check, ram_estimate);
stats.requested = true;

fprintf(['buildMobPeanutBigSparseStokes: assembling solve-grid ',...
    'sparse pair maps (mode=%s, P=%d, close pairs=%d, N_c=%d, ',...
    'N_check=%d, chunk pairs=%d).\n'],...
    build_mode, P, n_pairs, N_c, N_check, stats.chunk_pairs);

max_build_bytes = getOptField(opt,'big_sparse_max_build_bytes',inf);
if stats.big_sparse_build_bytes > max_build_bytes
    error('buildMobPeanutBigSparseStokes:BuildMemoryLimit',...
        ['Estimated sparse-entry build memory %.3g bytes exceeds ',...
         'opt.big_sparse_max_build_bytes %.3g bytes.'],...
        stats.big_sparse_build_bytes, max_build_bytes);
end

%% Early exit: no close pairs
if n_pairs == 0
    big_sparse              = emptyBigSparse(n_source_cols, n_u_rows,...
                                n_u_cols, pair_total_rows, n_rbm_rows, plan);
    big_sparse.ram_estimate = ram_estimate;
    stats.active            = true;
    stats.reason            = 'no_close_pairs';
    stats.build_time        = toc(timer);
    big_sparse.stats        = stats;
    return
end

%% Build sparse blocks
rout_base_c = rcheck(1:N_check) - q(1);
entries     = initSparseBuilders(P, N_c, N_check, n_pairs, opt, plan);

switch build_mode
    case 'precomputed'
        entries = appendPrecomputedBlocks(entries, geom, basis,...
                      rout_base_c, P_pair, plan);
    case 'streaming'
        entries = appendStreamingBlocks(entries, geom, basis,...
                      rout_base_c, P_pair, plan);
    otherwise
        error('buildMobPeanutBigSparseStokes:BadMode',...
            'Unsupported build mode "%s".', build_mode);
end

entries = flushAllSparseBuilders(entries);

%% Assemble output struct
big_sparse              = struct();
big_sparse.matrix_plan  = plan;
big_sparse.ram_estimate = ram_estimate;

if sparse_map_coarse
    big_sparse.M_source_corr = entries.source_corr.S;
else
    big_sparse.P_pair               = P_pair;
    big_sparse.source_scatter_rows  = buildSourceScatterRows(pairs, N_c, P);
    big_sparse.M_pair_nonp          = entries.pair_nonp.S;
end

big_sparse.M_rbm_corr = entries.rbm.S;

if direct_u_corr
    big_sparse.M_u_corr = entries.u.S;
else
    big_sparse.M_u_cross   = entries.u_cross.S;
    big_sparse.M_u_peanut  = entries.u_peanut.S;
    if sparse_map_coarse
        big_sparse.M_pair_proj = entries.pair_proj.S;
    end
end

%% Collect nnz stats
stats = collectNnzStats(stats, big_sparse, plan);
stats.active     = true;
stats.reason     = '';
stats.build_time = toc(timer);
big_sparse.stats = stats;

end % buildMobPeanutBigSparseStokes

% =========================================================================
%  BUILD-MODE & VALIDATION
% =========================================================================

function build_mode = resolveMobBigSparseBuildMode(opt)
build_mode = lower(char(getOptField(opt,'mob_big_sparse_build_mode',...
    'precomputed')));
if ~ismember(build_mode, {'precomputed','streaming'})
    error('buildMobPeanutBigSparseStokes:BadMode',...
        ['opt.mob_big_sparse_build_mode must be ''precomputed'' ',...
         'or ''streaming''.']);
end
end

function validateBigSparseInputs(geom, basis, N_check, build_mode)
opt = geom.opt;

% --- Required option flags ---
requiredFlags = {'cmap','self_correct','use_dense'};
for k = 1:numel(requiredFlags)
    if ~logical(getOptField(opt, requiredFlags{k}, false))
        error('buildMobPeanutBigSparseStokes:UnsupportedOption',...
            'opt.use_big_sparse=1 requires opt.%s=1.', requiredFlags{k});
    end
end

% --- Grid compatibility ---
if ~isequal(geom.rcheck, geom.rvec_out)
    error('buildMobPeanutBigSparseStokes:UnsupportedGrid',...
        'opt.use_big_sparse=1 is only supported on the solve grid.');
end
if N_check ~= round(N_check)
    error('buildMobPeanutBigSparseStokes:BadGridSize',...
        'rcheck length must be divisible by the number of particles.');
end

% --- Streaming-specific checks ---
if strcmp(build_mode,'streaming')
    if logical(getOptField(opt,'get_bndry_field',false))
        error('buildMobPeanutBigSparseStokes:StreamingBoundaryUnsupported',...
            ['opt.mob_big_sparse_build_mode=''streaming'' requires ',...
             'opt.get_bndry_field=0.']);
    end
    if logical(getOptField(opt,'reuse_pair_basis_by_sep',false))
        error('buildMobPeanutBigSparseStokes:StreamingPairReuseUnsupported',...
            ['opt.mob_big_sparse_build_mode=''streaming'' requires ',...
             'opt.reuse_pair_basis_by_sep=0.']);
    end
    if size(geom.pairs,1) > 0
        required_geom = {'rbase_in_f','rimage_vec','refine'};
        for k = 1:numel(required_geom)
            if ~isfield(geom,required_geom{k}) || isempty(geom.(required_geom{k}))
                error('buildMobPeanutBigSparseStokes:MissingStreamingGeom',...
                    'Streaming big sparse requires geom.%s.', required_geom{k});
            end
        end
        if ~isfield(basis,'Lc') || isempty(basis.Lc)
            error('buildMobPeanutBigSparseStokes:MissingLc',...
                'Streaming big sparse requires basis.Lc.');
        end
    end
end

% --- Precomputed-specific checks ---
if strcmp(build_mode,'precomputed') && size(geom.pairs,1) > 0
    if ~isfield(basis,'Cmap')    || isempty(basis.Cmap) ||...
       ~isfield(basis,'Cmap_FU') || isempty(basis.Cmap_FU)
        error('buildMobPeanutBigSparseStokes:MissingCmap',...
            'Precomputed big sparse requires actual per-pair Cmap and Cmap_FU data.');
    end
end

% --- Pair cache incompatibility ---
if size(geom.pairs,1) > 0 && isfield(basis,'pair_cache') &&...
        isfield(basis.pair_cache,'enabled') && basis.pair_cache.enabled
    error('buildMobPeanutBigSparseStokes:PairCacheUnsupported',...
        ['opt.use_big_sparse=1 currently requires actual per-pair maps. ',...
         'Set opt.reuse_pair_basis_by_sep=0 before building the pair basis.']);
end
end

% =========================================================================
%  STATS HELPERS
% =========================================================================

function stats = initBigSparseStats(n_pairs, N_c, N_check, ram_estimate)
counts = ram_estimate.counts;
plan   = ram_estimate.matrix_plan;
safe   = max(1, n_pairs);   % avoid divide-by-zero in per-pair ratios

stats = struct(...
    'requested',                false,...
    'active',                   false,...
    'backend',                  'global_block_sparse',...
    'matrix_plan',              plan,...
    'build_mode',               plan.build_mode,...
    'chunk_pairs',              plan.chunk_pairs,...
    'reason',                   'not_requested',...
    'n_pairs',                  n_pairs,...
    'N_c',                      N_c,...
    'N_check',                  N_check,...
    'used_pair_cache',          false,...
    'source_correction',        plan.source_correction,...
    'velocity_correction',      plan.velocity_correction,...
    'direct_u_corr',            plan.direct_u_corr,...
    'sparse_map_coarse',        plan.sparse_map_coarse,...
    'projector_mode',           plan.projector_mode,...
    'local_source_corr_entries',counts.source_corr / safe,...
    'local_pair_nonp_entries',  counts.pair_nonp   / safe,...
    'local_pair_proj_entries',  counts.pair_proj   / safe,...
    'local_u_entries',          counts.u           / safe,...
    'local_u_cross_entries',    counts.u_cross     / safe,...
    'local_u_peanut_entries',   counts.u_peanut    / safe,...
    'local_rbm_entries',        counts.rbm         / safe,...
    'nnz_source_corr',          0,...
    'nnz_u',                    0,...
    'nnz_u_cross',              0,...
    'nnz_u_peanut',             0,...
    'nnz_pair_nonp',            0,...
    'nnz_pair_proj',            0,...
    'nnz_rbm',                  0,...
    'nnz_source_scatter',       0,...
    'rotations_used',           false,...
    'big_sparse_matrix_bytes',  ram_estimate.big_sparse_matrix_bytes,...
    'big_sparse_auxiliary_bytes',ram_estimate.big_sparse_auxiliary_bytes,...
    'big_sparse_build_bytes',   ram_estimate.big_sparse_build_bytes,...
    'big_sparse_peak_bytes',    ram_estimate.big_sparse_peak_bytes,...
    'retained_pair_basis_bytes',ram_estimate.retained_pair_basis_bytes,...
    'solver_precompute_peak_bytes', ram_estimate.solver_precompute_peak_bytes,...
    'build_time',               0);
end

function stats = collectNnzStats(stats, big_sparse, plan)
%COLLECTNNZSTATS  Populate nnz fields from assembled big_sparse struct.
if plan.direct_u_corr
    stats.nnz_u = nnz(big_sparse.M_u_corr);
else
    stats.nnz_u_cross   = nnz(big_sparse.M_u_cross);
    stats.nnz_u_peanut  = nnz(big_sparse.M_u_peanut);
    if plan.sparse_map_coarse
        stats.nnz_pair_proj = nnz(big_sparse.M_pair_proj);
    end
end
if plan.sparse_map_coarse
    stats.nnz_source_corr = nnz(big_sparse.M_source_corr);
else
    stats.nnz_pair_nonp       = nnz(big_sparse.M_pair_nonp);
    stats.nnz_source_scatter  = numel(big_sparse.source_scatter_rows);
end
stats.nnz_rbm = nnz(big_sparse.M_rbm_corr);
end

% =========================================================================
%  EMPTY / ZERO-PAIR STRUCT
% =========================================================================

function big_sparse = emptyBigSparse(n_source_cols, n_u_rows, n_u_cols,...
    pair_total_rows, n_rbm_rows, plan)
big_sparse              = struct();
big_sparse.matrix_plan  = plan;

if plan.sparse_map_coarse
    big_sparse.M_source_corr = sparse(n_source_cols, n_source_cols);
else
    big_sparse.P_pair              = sparse(0,0);
    big_sparse.source_scatter_rows = zeros(0,1);
    big_sparse.M_pair_nonp         = sparse(pair_total_rows, n_source_cols);
end

big_sparse.M_rbm_corr = sparse(n_rbm_rows, n_source_cols);

if plan.direct_u_corr
    big_sparse.M_u_corr = sparse(n_u_rows, n_u_cols);
else
    big_sparse.M_u_cross  = sparse(n_u_rows, n_u_cols);
    big_sparse.M_u_peanut = sparse(n_u_rows, pair_total_rows);
    if plan.sparse_map_coarse
        big_sparse.M_pair_proj = sparse(pair_total_rows, n_source_cols);
    end
end
end

% =========================================================================
%  PAIR PROJECTION MATRIX
% =========================================================================

function P_pair = resolvePairProjectionMatrix(geom, basis)
opt = geom.opt;
if logical(getOptField(opt,'use_matrix_free_Lc_pair',true))
    P_pair = buildPairProjectionMatrix(geom.rbase_in_c(:));
    return
end

if isfield(basis,'Lc_pair') && ~isempty(basis.Lc_pair)
    P_pair = basis.Lc_pair;
elseif isfield(basis,'Lc') && ~isempty(basis.Lc)
    P_pair = getILpair(basis.Lc);
else
    error('buildMobPeanutBigSparseStokes:MissingDensePairProjector',...
        ['opt.use_matrix_free_Lc_pair=0 requires basis.Lc_pair ',...
         'or basis.Lc.']);
end

expected = 4 * geom.opt.N_c;
if ~isequal(size(P_pair), [expected expected])
    error('buildMobPeanutBigSparseStokes:BadDensePairProjector',...
        'basis.Lc_pair must have size %d-by-%d.', expected, expected);
end
end

function P_pair = buildPairProjectionMatrix(rbase)
K      = getKmat2D(rbase, 0);
L      = K * ((K'*K) \ K');
P_pair = getILpair(L);
end

% =========================================================================
%  SPARSE BUILDER INITIALISATION & FLUSHING
% =========================================================================

function entries = initSparseBuilders(P, N_c, N_check, n_pairs, opt, plan)
chunk_pairs = max(1, round(getOptField(opt,'mob_big_sparse_chunk_pairs',8)));
chunk_pairs = min(n_pairs, chunk_pairs);

n_coarse       = P * N_c;
n_source_cols  = 2 * n_coarse;
n_u_rows       = 2 * P * N_check;
n_u_cols       = 2 * n_coarse;
pair_rows      = 4 * N_c;
pair_total_rows = n_pairs * pair_rows;
n_rbm_rows     = 3 * P;

% Per-pair entry counts (used to pre-allocate builder buffers)
local_pair_nonp = pair_rows * (4*N_c);
local_u         = (4*N_check) * (4*N_c);
local_rbm       = 6 * (4*N_c);
chunk_nonp      = chunk_pairs * local_pair_nonp;
chunk_u         = chunk_pairs * local_u;
chunk_rbm       = chunk_pairs * local_rbm;

entries = struct();

% Source correction builders
if plan.sparse_map_coarse
    entries.source_corr = initSparseBuilder(chunk_nonp,  n_source_cols,  n_source_cols);
    entries.pair_nonp   = initSparseBuilder(0,            pair_total_rows, n_source_cols);
else
    entries.source_corr = initSparseBuilder(0,            n_source_cols,  n_source_cols);
    entries.pair_nonp   = initSparseBuilder(chunk_nonp,   pair_total_rows, n_source_cols);
end

entries.rbm = initSparseBuilder(chunk_rbm, n_rbm_rows, n_source_cols);

% Velocity correction builders
if plan.direct_u_corr
    entries.u         = initSparseBuilder(chunk_u, n_u_rows, n_u_cols);
    entries.u_cross   = initSparseBuilder(0,       n_u_rows, n_u_cols);
    entries.u_peanut  = initSparseBuilder(0,       n_u_rows, pair_total_rows);
    entries.pair_proj = initSparseBuilder(0,       pair_total_rows, n_source_cols);
else
    entries.u        = initSparseBuilder(0,       n_u_rows, n_u_cols);
    entries.u_cross  = initSparseBuilder(chunk_u, n_u_rows, n_u_cols);
    entries.u_peanut = initSparseBuilder(chunk_u, n_u_rows, pair_total_rows);
    if plan.sparse_map_coarse
        entries.pair_proj = initSparseBuilder(chunk_nonp, pair_total_rows, n_source_cols);
    else
        entries.pair_proj = initSparseBuilder(0,          pair_total_rows, n_source_cols);
    end
end
end

function entries = flushAllSparseBuilders(entries)
%FLUSHALLSPARSEBUILDERS Flush every builder field generically.
fields = fieldnames(entries);
for k = 1:numel(fields)
    entries.(fields{k}) = flushSparseBuilderBlock(entries.(fields{k}));
end
end

% =========================================================================
%  BLOCK ASSEMBLY — PRECOMPUTED & STREAMING
% =========================================================================

function entries = appendPrecomputedBlocks(entries, geom, basis,...
    rout_base_c, P_pair, plan)
%APPENDPRECOMPUTEDBLOCKS Assemble pair blocks from pre-built Cmap arrays.
pairs   = geom.pairs;
q       = geom.q(:);
N_c     = geom.opt.N_c;
N_check = numel(geom.rcheck) / numel(q);
P       = numel(q);
n_pairs = size(pairs,1);

for row = 1:n_pairs
    i = pairs(row,1);
    j = pairs(row,2);
    C_nonp   = basis.Cmap{i,j};
    Cmap_FU  = basis.Cmap_FU{i,j};
    [Ucross, Ecolloc] = buildStokesCoarsePairDense(...
        q, geom.rbase_in_c(:), rout_base_c, pairs, row);
    entries = appendPairBlocks(entries, pairs, row, N_c, N_check, P,...
        C_nonp, Cmap_FU, Ucross, Ecolloc, P_pair, plan);
end
end

function entries = appendStreamingBlocks(entries, geom, basis,...
    rout_base_c, P_pair, plan)
%APPENDSTREAMINGBLOCKS Assemble pair blocks by computing Cmap on the fly.
pairs       = geom.pairs;
q           = geom.q(:);
opt         = geom.opt;
N_c         = opt.N_c;
N_check     = numel(geom.rcheck) / numel(q);
P           = numel(q);
n_pairs     = size(pairs,1);
show_counter = logical(getOptField(opt,'show_counter',false));

svd_opts = struct(...
    'column_weight', logical(getOptField(opt,'column_weight',false)),...
    'left_weight',   logical(getOptField(opt,'left_weight',false)));

pair_opt               = opt;
pair_opt.project_force = true;
pair_opt.project       = true;
pair_opt.pair_basis_debug = false;

for row = 1:n_pairs
    pair = buildStokesMobilityPairData(q, geom.rbase_in_c(:),...
        geom.rbase_in_f(:), geom.rimage_vec, geom.refine, pairs,...
        pair_opt, basis.Lc, row, false, svd_opts, false, 'maps_only');
    [Ucross, Ecolloc] = buildStokesCoarsePairDense(...
        q, geom.rbase_in_c(:), rout_base_c, pairs, row);
    entries = appendPairBlocks(entries, pairs, row, N_c, N_check, P,...
        pair.Cmap, pair.Cmap_FU, Ucross, Ecolloc, P_pair, plan);
    if show_counter
        fprintf('buildMobPeanutBigSparseStokes: streamed pair %d/%d\n',...
            row, n_pairs);
    end
end
end

function entries = appendPairBlocks(entries, pairs, row, N_c, N_check, P,...
    C_nonp, Cmap_FU, Ucross, Ecolloc, P_pair, plan)
%APPENDPAIRBLOCKS Insert one pair's contribution into all sparse builders.
%
% plan.direct_u_corr and plan.sparse_map_coarse are resolved once from the
% mobility plan struct (passed in) rather than re-read from opt on every call.

i = pairs(row,1);
j = pairs(row,2);

in_idx       = pairCoarseInputIndices(i, j, N_c, P);
proj_out_idx = pairCoarseInputIndices(i, j, N_c, P);  % identical to in_idx layout
u_out_idx    = pairVelocityOutputIndices(i, j, N_check, P);
pair_idx     = (row-1)*(4*N_c)+1 : row*(4*N_c);
rbm_idx      = pairRigidOutputIndices(i, j, P);

% C_proj is needed for direct-u or sparse-map-coarse paths
if plan.sparse_map_coarse || plan.direct_u_corr
    C_proj = P_pair * C_nonp;
else
    C_proj = [];
end

% --- Source correction ---
if plan.sparse_map_coarse
    entries.source_corr = appendSparseBuilderBlock(entries.source_corr,...
        proj_out_idx, in_idx, C_proj);
else
    entries.pair_nonp = appendSparseBuilderBlock(entries.pair_nonp,...
        pair_idx, in_idx, C_nonp);
end

% --- Rigid-body motion correction ---
entries.rbm = appendSparseBuilderBlock(entries.rbm, rbm_idx, in_idx, -Cmap_FU);

% --- Velocity correction ---
if plan.direct_u_corr
    entries.u = appendSparseBuilderBlock(entries.u, u_out_idx, in_idx,...
        Ucross - Ecolloc*C_proj);
else
    entries.u_cross  = appendSparseBuilderBlock(entries.u_cross,...
        u_out_idx, in_idx, Ucross);
    entries.u_peanut = appendSparseBuilderBlock(entries.u_peanut,...
        u_out_idx, pair_idx, Ecolloc);
    if plan.sparse_map_coarse
        entries.pair_proj = appendSparseBuilderBlock(entries.pair_proj,...
            pair_idx, in_idx, C_proj);
    end
end
end

% =========================================================================
%  INDEX HELPERS
% =========================================================================

function rows = buildSourceScatterRows(pairs, N_c, P)
%BUILDSOURCESCATTERROWS Vectorised construction of scatter-row indices.
%
% For each pair (i,j) the 8*N_c destination rows are:
%   proj  x: body-i block, body-j block  (offset 0)
%   proj  y: body-i block, body-j block  (offset   n_coarse)
%   nonp  x: body-i block, body-j block  (offset 2*n_coarse)
%   nonp  y: body-i block, body-j block  (offset 3*n_coarse)

n_pairs  = size(pairs,1);
n_coarse = P * N_c;

block_offsets = [0; n_coarse; 2*n_coarse; 3*n_coarse];  % 4 x 1
local         = (0:N_c-1)';                               % N_c x 1

i_starts = (pairs(:,1)-1) * N_c;   % n_pairs x 1
j_starts = (pairs(:,2)-1) * N_c;   % n_pairs x 1

% ij_local: (2*N_c) x n_pairs  — body-i rows stacked above body-j rows
ij_local = [bsxfun(@plus, i_starts', local);    % N_c   x n_pairs
            bsxfun(@plus, j_starts', local)];   % 2*N_c x n_pairs

% Broadcast across 4 block offsets → (2*N_c) x n_pairs x 4, then flatten
% in column-major order: all rows for pair 1, then pair 2, …
rows = bsxfun(@plus, reshape(block_offsets,1,1,4), ij_local) + 1;
rows = reshape(permute(rows,[1 3 2]), 8*N_c*n_pairs, 1);
end

function idx = pairCoarseInputIndices(i, j, N_c, P)
%PAIRCOARSEINPUTINDICES Column indices for bodies i and j in [x;y] layout.
n_coarse = P * N_c;
idx_i    = (i-1)*N_c+1 : i*N_c;
idx_j    = (j-1)*N_c+1 : j*N_c;
idx      = [idx_i, idx_j, n_coarse+idx_i, n_coarse+idx_j]';
end

% pairProjectedSourceOutputIndices is identical to pairCoarseInputIndices;
% call through to avoid duplication.
function idx = pairProjectedSourceOutputIndices(i, j, N_c, P)
idx = pairCoarseInputIndices(i, j, N_c, P);
end

function idx = pairVelocityOutputIndices(i, j, N_check, P)
pm_check = P * N_check;
idx_i    = (i-1)*N_check+1 : i*N_check;
idx_j    = (j-1)*N_check+1 : j*N_check;
idx      = [idx_i, idx_j, pm_check+idx_i, pm_check+idx_j]';
end

function idx = pairRigidOutputIndices(i, j, P)
idx_i = (i-1)*3+1 : 3*i;
idx_j = (j-1)*3+1 : 3*j;
idx   = [idx_i, idx_j]';
if max(idx) > 3*P
    error('buildMobPeanutBigSparseStokes:BadRigidIndex',...
        'Internal rigid-motion output index out of range.');
end
end

% =========================================================================
%  SELF-TEST
% =========================================================================

function test_buildMobPeanutBigSparseStokes
fprintf(['buildMobPeanutBigSparseStokes self-test: random_discs_mc spy ',...
    'plots (and element magnitudes if plot_magn=1)\n']);

plot_magn = true;
close all;

if plot_magn
    P = 5;  N_c = 80;  N_f = 60;
    warning(['Check choices of N_c and N_f before inspecting matrix ',...
        'element magnitudes. Only small P is feasible with larger N_c.']);
else
    P = 20; N_c = 18;  N_f = 24;
end

opt           = get2Dparams(P, N_c, N_f);
opt.N_peanut  = 400 * plot_magn + 60 * ~plot_magn;

phi      = 0.65;
rad      = 1;
geom_opt = struct('domain','boxed','phi',phi,'rad',rad,...
    'min_gap',1e-3,'n_sweeps',200,'rng_seed',210421,'visualise',false);
[q, meta] = random_discs_mc(P, geom_opt);

opt.delta_pair              = 0.2;
opt.get_bndry_field         = 0;
opt.visualise_sol           = 0;
opt.visualise_grid          = 1;
opt.show_counter            = 0;
opt.cmap                    = 1;
opt.self_correct            = 1;
opt.use_dense               = 1;
opt.use_matrix_free_Lc_pair = 1;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp        = false;
opt.use_big_sparse          = true;

[geom, basis] = buildSelfTestMobilityData(q, opt);

% Four combinations of (mob_sparse_map_coarse, direct_u_corr)
cases = struct(...
    'sparse_map_coarse', {false, false, true,  true},...
    'direct_u_corr',     {true,  false, true,  false},...
    'label', {...
        'pair-based coarse rows, direct u',...
        'pair-based coarse rows, factored u',...
        'direct coarse, direct u  [SMALLEST storage & FASTEST apply]',...
        'direct coarse, factored u'});

fprintf('  P=%d, target phi=%.3f, actual phi=%.3f\n', P, phi, meta.phi);

for case_id = 1:numel(cases)
    geom_case = geom;
    geom_case.opt.mob_sparse_map_coarse   = cases(case_id).sparse_map_coarse;
    geom_case.opt.big_sparse_direct_u_corr = cases(case_id).direct_u_corr;

    [big_sparse, stats] = buildMobPeanutBigSparseStokes(geom_case, basis);

    fprintf(['  case %d: mob_sparse_map_coarse=%d direct_u_corr=%d ',...
        'pairs=%d  nnz: src_corr=%d pair_nonp=%d pair_proj=%d ',...
        'u=%d u_cross=%d u_peanut=%d rbm=%d\n'],...
        case_id, cases(case_id).sparse_map_coarse,...
        cases(case_id).direct_u_corr, stats.n_pairs,...
        stats.nnz_source_corr, stats.nnz_pair_nonp, stats.nnz_pair_proj,...
        stats.nnz_u, stats.nnz_u_cross, stats.nnz_u_peanut, stats.nnz_rbm);

    plotMobilityBigSparseSpyFigure(big_sparse, stats, P, meta.phi,...
        cases(case_id), plot_magn);
end

alignfigs;
end

% -------------------------------------------------------------------------

function [geom, basis] = buildSelfTestMobilityData(q, opt)
q   = q(:);
P   = numel(q);
rad = 1;
N_c  = opt.N_c;
N_f  = opt.N_f;
a_c  = opt.a_c;
Rp_c = opt.Rp_c;
Rp_f = opt.Rp_f;

% Create discretization
tout_c      = linspace(0, 2*pi, ceil(a_c*N_c)+1);
tout_c      = tout_c(1:end-1)';
rbase_out_c = rad * (cos(tout_c) + 1i*sin(tout_c));

tin_c      = linspace(0, 2*pi, N_c+1)'; tin_c = tin_c(1:end-1);
tin_f      = linspace(0, 2*pi, N_f+1)'; tin_f = tin_f(1:end-1);
rbase_in_c = Rp_c * exp(1i*tin_c);
rbase_in_f = Rp_f * exp(1i*tin_f);

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q, opt);

% Assemble global coarse input and output grids
rvec_in_c = zeros(P*N_c, 1);
rout      = zeros(P*numel(rbase_out_c), 1);
for k = 1:P
    rvec_in_c((k-1)*N_c+1:k*N_c)                          = q(k) + rbase_in_c;
    rout((k-1)*numel(rbase_out_c)+1:k*numel(rbase_out_c)) = q(k) + rbase_out_c;
end

svd_opts = struct('column_weight',false,'left_weight',false);
[U, Y, Lc] = getSelfPseudoMobilityStokes(...
    1, q, rbase_in_c, rbase_out_c, [], [0, ceil(a_c*N_c)], svd_opts);

pair_opt                  = opt;
pair_opt.project_force    = true;
pair_opt.project          = true;
pair_opt.pair_basis_debug = 0;
pair_opt.rad              = ones(P,1);

[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU,pair_cache] =...
    getPairBasisStokes(q, rbase_in_c, rbase_in_f, rimage_vec, refine,...
    pairs, pair_opt, Lc{1}, rbase_out_c, svd_opts);

% --- geom ---
geom               = struct();
geom.rbase_in_c    = rbase_in_c;
geom.rbase_in_f    = rbase_in_f;
geom.rvec_in       = rvec_in_c;
geom.rimage_vec    = rimage_vec;
geom.opt           = pair_opt;       
geom.opt.get_bndry_field  = 0;
geom.opt.parallel_solve   = false;
geom.rvec_out      = rout;
geom.rcheck        = rout;
geom.q             = q;
geom.pairs         = pairs;
geom.refine        = refine;

% --- basis ---
basis              = struct();
basis.U            = U;
basis.Y            = Y;
basis.Lc           = Lc{1};
basis.Upf          = UB_all;
basis.Ypf          = YB_all;
basis.DC_all       = UC_all;
basis.YC_all       = YC_all;
basis.Cmap         = Cmap;
basis.Cmap_FU      = Cmap_FU;
basis.Lc_pair      = getILpair(Lc{1});
basis.pair_cache   = pair_cache;
end

% -------------------------------------------------------------------------

function plotMobilityBigSparseSpyFigure(big_sparse, stats, P, phi,...
    case_info, plot_magn)

% Collect matrices and titles based on active mode
matrices = {};
titles   = {};

if case_info.sparse_map_coarse
    matrices{end+1} = big_sparse.M_source_corr;
    titles{end+1}   = 'Source correction map (direct coarse)';
else
    matrices{end+1} = big_sparse.M_pair_nonp;
    titles{end+1}   = 'Source correction map (pair non-projected)';
end

if case_info.direct_u_corr
    figure_name     = 'buildMobPeanutBigSparseStokes: direct maps';
    matrices{end+1} = big_sparse.M_u_corr;
    titles{end+1}   = 'Velocity correction map (direct)';
else
    figure_name = 'buildMobPeanutBigSparseStokes: factored maps';
    if case_info.sparse_map_coarse
        matrices{end+1} = big_sparse.M_pair_proj;
        titles{end+1}   = 'Projected pair map (for factored u)';
    end
    matrices{end+1} = big_sparse.M_u_cross;
    titles{end+1}   = 'Velocity correction: cross part';
    matrices{end+1} = big_sparse.M_u_peanut;
    titles{end+1}   = 'Velocity correction: peanut part';
end

[nrows, ncols] = chooseSpyTileShape(numel(matrices));
sup_title = {...
    sprintf('P=%d, phi=%.2f, pairs=%d, N_c=%d', P, phi, stats.n_pairs, stats.N_c),...
    case_info.label,...
    '(RBM correction not shown)'};

figure('Name', sprintf('%s (%s)', figure_name, case_info.label), 'Color','w');
tiledlayout(nrows, ncols, 'TileSpacing','compact', 'Padding','compact');

for k = 1:numel(matrices)
    nexttile;
    if plot_magn
        imagesc(log10(abs(full(matrices{k}))));
        colorbar;
        axis equal tight;
        title([titles{k}, ' [log_{10}|·|]'], 'Interpreter','tex');
    else
        spy(matrices{k});
        title(titles{k}, 'Interpreter','none');
    end
end

if plot_magn
    sup_title{end+1} =...
        'WARNING: large magnitudes may indicate small N_c';
end
sgtitle(sup_title);
end

function [nrows, ncols] = chooseSpyTileShape(nplots)
nrows = 1 + (nplots > 3);
ncols = ceil(nplots / nrows);
end
