function [big_sparse,stats] = buildResPeanutBigSparseStokes(geom,basis)
%BUILDRESPEANUTBIGSPARSESTOKES Build solve-grid resistance pair maps. These 
% are both to correct coarse sources and to correct the velocity 
% field on the close pairs themselves due to coarse sources.
%
% Calling with no inputs runs a small visual self-test and shows the sparse
% matrix structures with spy and, if requested, the actual magnitudes of
% the compressed matrices.
%
% opt.res_sparse_map_coarse, opt.res_big_sparse_u_corr_mode, and
% opt.res_smaller_ucorr select the matrix set built here:
%   res_sparse_map_coarse = 0:
%       Coarse source corrections use pair-based sorting of coarse sources.
%   res_sparse_map_coarse = 1:
%       Coarse source corrections directly map projected one-body coarse
%       sources to projected coarse-source increments (per particle, not
%       per pair).
%   res_big_sparse_u_corr_mode = 'combined':
%       res_smaller_ucorr = 0  builds a combined velocity correction map
%         acting on both self and pair coarse sources.
%       res_smaller_ucorr = 1  builds a structured dense payload for
%         pair-dependent cross corrections, reusing the self block during
%         application. The legacy alias opt.smaller_ucorr is also accepted.
%   res_big_sparse_u_corr_mode = 'factored':
%       Velocity corrections are factored into cross and pair components.
%
% % Note: The fastest and least memory-heavy setup is
% res_big_sparse_u_corr_mode = 'combined', with res_small_ucorr = 1 and 
% res_sparse_map_coarse = 0. Setting big_sparse_u_corr_mode = 'factored' might be
% slighly better behaved numerically, but the difference seems to be tiny.
%
% All the options described here are compared in the script
% apr24_resistance_sparse_map_coarse_compare.m in the experiments/ folder.

if nargin == 0
    test_buildResPeanutBigSparseStokes;
    big_sparse = [];
    stats      = [];
    return
end

timer = tic;

opt     = geom.opt;
q       = geom.q(:);
pairs   = geom.pairs;
rcheck  = geom.rcheck;
P       = opt.P;
N_c     = opt.N_c;
N_check = numel(rcheck) / P;
n_pairs = size(pairs,1);
n_coarse = P*N_c;

%% Dimension constants (derived once, shared throughout)
pair_rows       = 4*N_c;
pair_total_rows = n_pairs * pair_rows;
n_source_cols   = 2*n_coarse;
n_global_u_rows = 2*P*N_check;

validateBigSparseInputs(geom, basis, N_check);

%% RAM plan & stats
ram_estimate = estimateResPeanutBigSparseRamStokes(P, N_c, N_check, n_pairs, opt);
plan         = ram_estimate.matrix_plan;
stats        = initBigSparseStats(n_pairs, N_c, N_check, ram_estimate);
stats.requested = true;

fprintf(['buildResPeanutBigSparseStokes: assembling solve-grid ',...
    'sparse pair maps (P=%d, close pairs=%d, N_c=%d, ',...
    'N_check=%d, velocity=%s, chunk pairs=%d).\n'],...
    P, n_pairs, N_c, N_check, stats.velocity_correction, stats.chunk_pairs);

%% Build-memory guard
max_build_bytes = getOptField(opt,'res_big_sparse_max_build_bytes',...
    getOptField(opt,'big_sparse_max_build_bytes',inf));
if stats.big_sparse_build_bytes > max_build_bytes
    error('buildResPeanutBigSparseStokes:BuildMemoryLimit',...
        ['Estimated sparse-entry build memory %.3g bytes exceeds ',...
         'the configured limit %.3g bytes.'],...
        stats.big_sparse_build_bytes, max_build_bytes);
end

%% Early exit: no close pairs
if n_pairs == 0
    big_sparse = emptyBigSparse(pair_total_rows, n_source_cols,...
        n_global_u_rows, N_c, N_check, ram_estimate);
    stats.active     = true;
    stats.reason     = 'no_close_pairs';
    stats.build_time = toc(timer);
    big_sparse.stats = stats;
    return
end

%% Assemble one dense pair contribution at a time into global sparse maps
rout_base_c = rcheck(1:N_check) - q(1);
entries     = initSparseBuilders(P, N_c, N_check, n_pairs, opt, plan);

% Cross-velocity map is not needed for the combined+smaller path
need_cross_velocity_map = ~(strcmp(plan.u_corr_mode,'combined') &&...
    plan.smaller_u_corr);

for row = 1:n_pairs
    i = pairs(row,1);
    j = pairs(row,2);
    C_pair  = basis.Cmap{i,j};
    Cmap_FU = basis.Cmap_FU{i,j};

    if need_cross_velocity_map
        [Ucross, Ecolloc] = buildStokesCoarsePairDense(...
            q, geom.rbase_in_c(:), rout_base_c, pairs, row);
    else
        Ucross  = [];
        Ecolloc = buildStokesPairEcolloc(...
            q, geom.rbase_in_c(:), rout_base_c, pairs, row);
    end

    entries = appendPairBlocks(entries, pairs, row, N_c, N_check, P,...
        C_pair, Cmap_FU, Ucross, Ecolloc, plan);
end

entries = flushAllSparseBuilders(entries);

%% Assemble output struct
big_sparse              = struct();
big_sparse.matrix_plan  = plan;
big_sparse.ram_estimate = ram_estimate;
big_sparse.N_c          = N_c;
big_sparse.pairs        = pairs;

if plan.sparse_map_coarse
    big_sparse.M_source_corr = entries.source_corr.S;
else
    big_sparse.source_scatter_rows = buildSourceScatterRows(pairs, N_c, P);
end

big_sparse.velocity_scatter_rows = buildVelocityScatterRows(pairs, N_check, P);
big_sparse.ft_scatter_rows       = buildForceTorqueScatterRows(pairs, P);
big_sparse.M_pair                = entries.pair.S;

switch plan.u_corr_mode
    case 'combined'
        if plan.smaller_u_corr
            big_sparse.u_corr_small = entries.u_small;
        else
            big_sparse.M_u_corr = entries.u.S;
        end
    otherwise
        big_sparse.M_u_cross  = entries.u_cross.S;
        big_sparse.M_u_peanut = entries.u_peanut.S;
end

big_sparse.M_ft_corr = entries.ft.S;

%% Collect nnz stats
stats            = collectNnzStats(stats, big_sparse, plan);
stats.active     = true;
stats.reason     = '';
stats.build_time = toc(timer);
big_sparse.stats = stats;

end % buildResPeanutBigSparseStokes

% =========================================================================
%  VALIDATION
% =========================================================================

function validateBigSparseInputs(geom, basis, N_check)
opt = geom.opt;

% Required option flags
if ~logical(getOptField(opt,'cmap',false))
    error('buildResPeanutBigSparseStokes:UnsupportedOption',...
        'opt.use_big_sparse=1 requires opt.cmap=1.');
end
if ~logical(getOptField(opt,'self_correct',false))
    error('buildResPeanutBigSparseStokes:UnsupportedOption',...
        'opt.use_big_sparse=1 requires opt.self_correct=1.');
end
if getOptField(opt,'N_peanut',0) <= 0
    error('buildResPeanutBigSparseStokes:UnsupportedOption',...
        'opt.use_big_sparse=1 requires opt.N_peanut > 0.');
end

% Grid compatibility
if ~isequal(geom.rcheck, geom.rvec_out)
    error('buildResPeanutBigSparseStokes:UnsupportedGrid',...
        'Resistance big sparse pair maps are only built on the solve grid.');
end
if N_check ~= round(N_check)
    error('buildResPeanutBigSparseStokes:BadGridSize',...
        'rcheck length must be divisible by the number of particles.');
end

% Per-pair data requirements
if size(geom.pairs,1) > 0
    if ~isfield(basis,'Cmap')    || isempty(basis.Cmap) ||...
       ~isfield(basis,'Cmap_FU') || isempty(basis.Cmap_FU)
        error('buildResPeanutBigSparseStokes:MissingPairMaps',...
            ['Resistance big sparse assembly requires actual per-pair ',...
             'Cmap and Cmap_FU data.']);
    end
    if isfield(basis,'pair_cache') && isfield(basis.pair_cache,'enabled')...
            && basis.pair_cache.enabled
        error('buildResPeanutBigSparseStokes:PairCacheUnsupported',...
            ['opt.use_big_sparse=1 requires actual per-pair maps. ',...
             'Set opt.reuse_pair_basis_by_sep=0 before building the pair basis.']);
    end
end
end

% =========================================================================
%  STATS HELPERS
% =========================================================================

function stats = initBigSparseStats(n_pairs, N_c, N_check, ram_estimate)
counts = ram_estimate.counts;
plan   = ram_estimate.matrix_plan;
safe   = max(1, n_pairs);

stats = struct(...
    'requested',                false,...
    'active',                   false,...
    'backend',                  'global_block_sparse',...
    'reason',                   'not_requested',...
    'matrix_plan',              plan,...
    'chunk_pairs',              plan.chunk_pairs,...
    'n_pairs',                  n_pairs,...
    'N_c',                      N_c,...
    'N_check',                  N_check,...
    'used_pair_cache',          false,...
    'source_correction',        plan.source_correction,...
    'sparse_map_coarse',        plan.sparse_map_coarse,...
    'velocity_correction',      plan.velocity_correction,...
    'u_corr_mode',              plan.u_corr_mode,...
    'combined_u_corr',          plan.combined_u_corr,...
    'smaller_u_corr',           plan.smaller_u_corr,...
    'force_torque_correction',  plan.force_torque_correction,...
    'local_source_corr_entries',counts.source_corr / safe,...
    'local_pair_entries',       counts.pair        / safe,...
    'local_u_entries',          counts.u           / safe,...
    'local_u_small_entries',    counts.u_small     / safe,...
    'local_u_cross_entries',    counts.u_cross     / safe,...
    'local_u_peanut_entries',   counts.u_peanut    / safe,...
    'local_ft_entries',         counts.ft          / safe,...
    'nnz_source_corr',          0,...
    'nnz_pair',                 0,...
    'nnz_u',                    0,...
    'nnz_u_small',              0,...
    'nnz_u_cross',              0,...
    'nnz_u_peanut',             0,...
    'nnz_ft',                   0,...
    'nnz_source_scatter',       0,...
    'nnz_velocity_scatter',     0,...
    'nnz_ft_scatter',           0,...
    'big_sparse_matrix_bytes',  ram_estimate.big_sparse_matrix_bytes,...
    'big_sparse_auxiliary_bytes',ram_estimate.big_sparse_auxiliary_bytes,...
    'big_sparse_build_bytes',   ram_estimate.big_sparse_build_bytes,...
    'big_sparse_peak_bytes',    ram_estimate.big_sparse_peak_bytes,...
    'build_time',               0);
end

function stats = collectNnzStats(stats, big_sparse, plan)
%COLLECTNNZSTATS Populate nnz fields from the assembled big_sparse struct.
stats.nnz_pair = nnz(big_sparse.M_pair);
stats.nnz_ft   = nnz(big_sparse.M_ft_corr);

switch plan.u_corr_mode
    case 'combined'
        if plan.smaller_u_corr
            stats.nnz_u_small = numel(big_sparse.u_corr_small.cross_ij) +...
                                 numel(big_sparse.u_corr_small.cross_ji);
        else
            stats.nnz_u = nnz(big_sparse.M_u_corr);
        end
    otherwise
        stats.nnz_u_cross  = nnz(big_sparse.M_u_cross);
        stats.nnz_u_peanut = nnz(big_sparse.M_u_peanut);
end

if plan.sparse_map_coarse
    stats.nnz_source_corr = nnz(big_sparse.M_source_corr);
else
    stats.nnz_source_scatter = numel(big_sparse.source_scatter_rows);
end

stats.nnz_velocity_scatter = numel(big_sparse.velocity_scatter_rows);
stats.nnz_ft_scatter       = numel(big_sparse.ft_scatter_rows);
end

% =========================================================================
%  EMPTY / ZERO-PAIR STRUCT
% =========================================================================

function big_sparse = emptyBigSparse(pair_total_rows, n_source_cols,...
    n_global_u_rows, N_c, N_check, ram_estimate)
plan = ram_estimate.matrix_plan;

big_sparse               = struct();
big_sparse.matrix_plan   = plan;
big_sparse.ram_estimate  = ram_estimate;
big_sparse.N_c           = [];
big_sparse.pairs         = zeros(0,2);

if plan.sparse_map_coarse
    big_sparse.M_source_corr = sparse(n_source_cols, n_source_cols);
else
    big_sparse.source_scatter_rows = zeros(0,1);
end

big_sparse.velocity_scatter_rows = zeros(0,1);
big_sparse.ft_scatter_rows       = zeros(0,1);
big_sparse.M_pair                = sparse(pair_total_rows, n_source_cols);

switch plan.u_corr_mode
    case 'combined'
        if plan.smaller_u_corr
            big_sparse.u_corr_small = initSmallUCorrBuilder(...
                N_c, N_check, 0, plan.chunk_pairs);
        else
            big_sparse.M_u_corr = sparse(0, n_source_cols + pair_total_rows);
        end
    otherwise
        big_sparse.M_u_cross  = sparse(0, n_source_cols);
        big_sparse.M_u_peanut = sparse(0, pair_total_rows);
end

big_sparse.M_ft_corr       = sparse(0, n_source_cols);
big_sparse.n_global_u_rows = n_global_u_rows;
end

% =========================================================================
%  SPARSE BUILDER INITIALISATION & FLUSHING
% =========================================================================

function entries = initSparseBuilders(P, N_c, N_check, n_pairs, opt, plan)
chunk_pairs = max(1, round(getOptField(opt,'res_big_sparse_chunk_pairs',8)));
chunk_pairs = min(n_pairs, chunk_pairs);

n_coarse        = P * N_c;
n_source_cols   = 2 * n_coarse;
pair_rows       = 4 * N_c;
pair_total_rows = n_pairs * pair_rows;
pair_u_rows     = n_pairs * 4 * N_check;
pair_ft_rows    = n_pairs * 6;

% Per-pair entry counts for buffer pre-allocation
local_pair = pair_rows * (4*N_c);
local_u    = (4*N_check) * (4*N_c);
local_ft   = 6 * (4*N_c);
chunk_pair = chunk_pairs * local_pair;
chunk_u    = chunk_pairs * local_u;
chunk_ft   = chunk_pairs * local_ft;

u_corr_mode   = plan.u_corr_mode;
smaller_u_corr = plan.smaller_u_corr;

entries = struct();

% Source correction
if plan.sparse_map_coarse
    entries.source_corr = initSparseBuilder(chunk_pair, n_source_cols, n_source_cols);
else
    entries.source_corr = initSparseBuilder(0,          n_source_cols, n_source_cols);
end

entries.pair = initSparseBuilder(chunk_pair, pair_total_rows, n_source_cols);

% Velocity correction
switch u_corr_mode
    case 'combined'
        if smaller_u_corr
            entries.u       = initSparseBuilder(0, pair_u_rows, n_source_cols);
            entries.u_small = initSmallUCorrBuilder(N_c, N_check, n_pairs, chunk_pairs);
        else
            entries.u       = initSparseBuilder(chunk_pairs*2*local_u, pair_u_rows,...
                                  n_source_cols + pair_total_rows);
            entries.u_small = initSmallUCorrBuilder(N_c, N_check, 0, chunk_pairs);
        end
        entries.u_cross  = initSparseBuilder(0, pair_u_rows, n_source_cols);
        entries.u_peanut = initSparseBuilder(0, pair_u_rows, pair_total_rows);
    otherwise
        entries.u        = initSparseBuilder(0,       pair_u_rows, n_source_cols);
        entries.u_small  = initSmallUCorrBuilder(N_c, N_check, 0, chunk_pairs);
        entries.u_cross  = initSparseBuilder(chunk_u, pair_u_rows, n_source_cols);
        entries.u_peanut = initSparseBuilder(chunk_u, pair_u_rows, pair_total_rows);
end

entries.ft = initSparseBuilder(chunk_ft, pair_ft_rows, n_source_cols);
end

function entries = flushAllSparseBuilders(entries)
%FLUSHALLSPARSEBUILDERS Flush all sparse builder fields.
% u_small uses a different (dense) flush; all others use flushSparseBuilderBlock.
sparse_fields = {'source_corr','pair','u','u_cross','u_peanut','ft'};
for k = 1:numel(sparse_fields)
    f = sparse_fields{k};
    entries.(f) = flushSparseBuilderBlock(entries.(f));
end
entries.u_small = flushSmallUCorrBuilder(entries.u_small);
end

% =========================================================================
%  PAIR BLOCK ASSEMBLY
% =========================================================================

function entries = appendPairBlocks(entries, pairs, row, N_c, N_check, P,...
    C_pair, Cmap_FU, Ucross, Ecolloc, plan)
%APPENDPAIRBLOCKS Insert one pair's contribution into all sparse builders.
%
% plan is passed in (resolved once per build) so no per-pair opt lookups
% are needed here.

i = pairs(row,1);
j = pairs(row,2);

in_idx     = pairCoarseInputIndices(i, j, N_c, P);
source_idx = in_idx;                               % identical layout
pair_idx   = (row-1)*(4*N_c)+1 : row*(4*N_c);
u_idx      = pairVelocityBlockIndices(row, N_check);
ft_idx     = pairForceTorqueBlockIndices(row);

% Source correction
if plan.sparse_map_coarse
    entries.source_corr = appendSparseBuilderBlock(entries.source_corr,...
        source_idx, in_idx, C_pair);
end

entries.pair = appendSparseBuilderBlock(entries.pair, pair_idx, in_idx, C_pair);

% Velocity correction
switch plan.u_corr_mode
    case 'combined'
        if plan.smaller_u_corr
            entries.u_small = appendSmallUCorrBlocks(entries.u_small,...
                row, N_c, N_check, Ecolloc);
        else
            n_source_cols = 2 * P * N_c;
            entries.u = appendSparseBuilderBlock(entries.u, u_idx, in_idx, Ucross);
            entries.u = appendSparseBuilderBlock(entries.u, u_idx,...
                n_source_cols + pair_idx, -Ecolloc);
        end
    otherwise
        entries.u_cross  = appendSparseBuilderBlock(entries.u_cross,...
            u_idx, in_idx, Ucross);
        entries.u_peanut = appendSparseBuilderBlock(entries.u_peanut,...
            u_idx, pair_idx, Ecolloc);
end

% Force/torque correction
entries.ft = appendSparseBuilderBlock(entries.ft, ft_idx, in_idx, Cmap_FU);
end

% =========================================================================
%  SMALL U-CORR DENSE PAYLOAD
% =========================================================================

function payload = initSmallUCorrBuilder(N_c, N_check, n_pairs, chunk_pairs)
payload = struct(...
    'cross_ij',   zeros(2*N_check, 2*N_c, n_pairs),...
    'cross_ji',   zeros(2*N_check, 2*N_c, n_pairs),...
    'chunk_pairs', chunk_pairs);
end

function payload = appendSmallUCorrBlocks(payload, row, N_c, N_check, Ecolloc)
[payload.cross_ij(:,:,row), payload.cross_ji(:,:,row)] =...
    extractPairCrossBlocks(Ecolloc, N_c, N_check);
end

function payload = flushSmallUCorrBuilder(payload)
% Dense small-u_corr payload is assembled directly into its final arrays;
% nothing to flush.
end

function [cross_ij, cross_ji] = extractPairCrossBlocks(Ecolloc, N_c, N_check)
%EXTRACTPAIRCROSSBLOCKS Extract the off-diagonal body blocks from Ecolloc.
%
% Ecolloc is (4*N_check) x (4*N_c), laid out as:
%   rows: [u_i_x | u_j_x | u_i_y | u_j_y]  (each block of N_check rows)
%   cols: [f_i_x | f_j_x | f_i_y | f_j_y]  (each block of N_c cols)
%
% cross_ij = effect of body-j sources on body-i velocities
% cross_ji = effect of body-i sources on body-j velocities

rows_i = [1:N_check,   2*N_check+1:3*N_check];   % x and y rows for body i
rows_j = [N_check+1:2*N_check, 3*N_check+1:4*N_check];
cols_i = [1:N_c,       2*N_c+1:3*N_c];           % x and y cols for body i
cols_j = [N_c+1:2*N_c, 3*N_c+1:4*N_c];

cross_ij = Ecolloc(rows_i, cols_j);
cross_ji = Ecolloc(rows_j, cols_i);
end

% =========================================================================
%  ECOLLOC-ONLY BUILD (smaller_u_corr path)
% =========================================================================

function Ecolloc = buildStokesPairEcolloc(q, rbase_in_c, rout_base_c, pairs, row)
i = pairs(row,1);
j = pairs(row,2);
rin_pair_c = [q(i)+rbase_in_c(:); q(j)+rbase_in_c(:)];
rout_pair  = [q(i)+rout_base_c(:); q(j)+rout_base_c(:)];
Ecolloc    = stokSLPmat(rin_pair_c, rout_pair, 1);
end

% =========================================================================
%  SCATTER ROW BUILDERS — vectorised
% =========================================================================

function rows = buildSourceScatterRows(pairs, N_c, P)
%BUILDSOURCESCATTERROWS Vectorised: source scatter indices for all pairs.
% Each pair contributes 4*N_c rows (x and y blocks for bodies i and j).
n_pairs  = size(pairs,1);
n_coarse = P * N_c;

% Block offsets: [x-block, y-block] = [0, n_coarse]
block_offsets = [0; n_coarse];          % 2 x 1
local         = (0:N_c-1)';            % N_c x 1

i_starts = (pairs(:,1)-1) * N_c;       % n_pairs x 1
j_starts = (pairs(:,2)-1) * N_c;

% ij_local: (2*N_c) x n_pairs
ij_local = [bsxfun(@plus, i_starts', local);
            bsxfun(@plus, j_starts', local)];

% Broadcast across 2 block offsets → (2*N_c) x n_pairs x 2, then flatten
rows = bsxfun(@plus, reshape(block_offsets,1,1,2), ij_local) + 1;
rows = reshape(permute(rows,[1 3 2]), 4*N_c*n_pairs, 1);
end

function rows = buildVelocityScatterRows(pairs, N_check, P)
%BUILDVELOCITYSCATTERROWS Vectorised: velocity scatter indices for all pairs.
% Each pair contributes 4*N_check rows.
n_pairs  = size(pairs,1);
pm_check = P * N_check;

block_offsets = [0; pm_check];
local         = (0:N_check-1)';

i_starts = (pairs(:,1)-1) * N_check;
j_starts = (pairs(:,2)-1) * N_check;

ij_local = [bsxfun(@plus, i_starts', local);
            bsxfun(@plus, j_starts', local)];

rows = bsxfun(@plus, reshape(block_offsets,1,1,2), ij_local) + 1;
rows = reshape(permute(rows,[1 3 2]), 4*N_check*n_pairs, 1);
end

function rows = buildForceTorqueScatterRows(pairs, P)
%BUILDFORCETORQUESCATTERROWS Vectorised: FT scatter indices for all pairs.
% Each pair contributes 6 rows (3 DOF x 2 bodies).
n_pairs = size(pairs,1);
local   = (0:2)';                       % 3 x 1

i_starts = (pairs(:,1)-1) * 3;         % n_pairs x 1
j_starts = (pairs(:,2)-1) * 3;

% ij_local: 6 x n_pairs
ij_local = [bsxfun(@plus, i_starts', local);
            bsxfun(@plus, j_starts', local)];

rows = reshape(ij_local + 1, 6*n_pairs, 1);

if max(rows) > 3*P
    error('buildResPeanutBigSparseStokes:BadForceTorqueIndex',...
        'Internal force/torque scatter index out of range.');
end
end

% =========================================================================
%  INDEX HELPERS
% =========================================================================

function idx = pairCoarseInputIndices(i, j, N_c, P)
%PAIRCOARSEINPUTINDICES Column indices for bodies i and j in [x;y] layout.
% Also used as source output indices (layouts are identical).
n_coarse = P * N_c;
idx_i    = (i-1)*N_c+1 : i*N_c;
idx_j    = (j-1)*N_c+1 : j*N_c;
idx      = [idx_i, idx_j, n_coarse+idx_i, n_coarse+idx_j]';
end

% pairSourceOutputIndices has the same layout as pairCoarseInputIndices.
% It is aliased here for call-site clarity.
function idx = pairSourceOutputIndices(i, j, N_c, P)
idx = pairCoarseInputIndices(i, j, N_c, P);
end

function idx = pairVelocityOutputIndices(i, j, N_check, P)
pm_check = P * N_check;
idx_i    = (i-1)*N_check+1 : i*N_check;
idx_j    = (j-1)*N_check+1 : j*N_check;
idx      = [idx_i, idx_j, pm_check+idx_i, pm_check+idx_j]';
end

function idx = pairForceTorqueOutputIndices(i, j, P)
idx_i = (i-1)*3+1 : 3*i;
idx_j = (j-1)*3+1 : 3*j;
idx   = [idx_i, idx_j]';
if max(idx) > 3*P
    error('buildResPeanutBigSparseStokes:BadForceTorqueIndex',...
        'Internal force/torque output index out of range.');
end
end

function idx = pairVelocityBlockIndices(row, N_check)
idx = ((row-1)*(4*N_check)+1 : row*(4*N_check))';
end

function idx = pairForceTorqueBlockIndices(row)
idx = ((row-1)*6+1 : row*6)';
end

% =========================================================================
%  SELF-TEST
% =========================================================================

function test_buildResPeanutBigSparseStokes
fprintf(['buildResPeanutBigSparseStokes self-test: spy plots ',...
    '(and imagesc for element magnitudes if plot_magn=1)\n']);

plot_magn = true;
close all;

if plot_magn
    P = 5;  N_c = 100;  N_f = 60;
    warning(['Check N_c and N_f before inspecting matrix element magnitudes. ',...
        'Only small P is feasible with larger N_c due to RAM constraints.']);
else
    P = 20; N_c = 18;   N_f = 24;
end

opt          = get2Dparams(P, N_c, N_f);
opt.N_peanut = 400 * plot_magn + 60 * ~plot_magn;

phi      = 0.65;
rad      = 1;
geom_opt = struct('domain','boxed','phi',phi,'rad',rad,...
    'min_gap',1e-3,'n_sweeps',200,'rng_seed',210421,'visualise',false);
[q, meta] = random_discs_mc(P, geom_opt);

opt.delta_pair              = 0.2;
opt.get_bndry_field         = 0;
opt.visualise_grid          = 1;
opt.show_counter            = 1;
opt.cmap                    = 1;
opt.self_correct            = 1;
opt.use_dense               = 1;
opt.reuse_pair_basis_by_sep = false;
opt.use_fmm                 = false;

[geom, basis] = buildSelfTestData(q, opt);

cases = struct(...
    'sparse_map_coarse', {false, false, true,  true},...
    'u_corr_mode',       {'combined','factored','combined','factored'},...
    'label', {...
        'pair-based coarse sorting, combined u',...
        'pair-based coarse sorting, factored u',...
        'direct coarse, combined u',...
        'direct coarse, factored u'});

fprintf('  P=%d, target phi=%.3f, actual phi=%.3f\n', P, phi, meta.phi);

for case_id = 3:numel(cases)
    geom_case = geom;
    geom_case.opt.res_sparse_map_coarse      = cases(case_id).sparse_map_coarse;
    geom_case.opt.res_big_sparse_u_corr_mode = cases(case_id).u_corr_mode;

    [big_sparse, stats] = buildResPeanutBigSparseStokes(geom_case, basis);

    fprintf(['  case %d: res_sparse_map_coarse=%d u_mode=%s pairs=%d ',...
        'nnz: src_corr=%d pair=%d u=%d u_cross=%d u_peanut=%d ft=%d\n'],...
        case_id, cases(case_id).sparse_map_coarse, cases(case_id).u_corr_mode,...
        stats.n_pairs, stats.nnz_source_corr, stats.nnz_pair,...
        stats.nnz_u, stats.nnz_u_cross, stats.nnz_u_peanut, stats.nnz_ft);

    plotResBigSparseSpyFigure(big_sparse, stats, P, meta.phi,...
        cases(case_id), plot_magn);
end

alignfigs;
end

% -------------------------------------------------------------------------

function [geom, basis] = buildSelfTestData(q, opt)
q   = q(:);
P   = numel(q);
N_c = opt.N_c;
N_f = opt.N_f;
a_c = opt.a_c;

% Discretize
tout_c      = linspace(0, 2*pi, ceil(a_c*N_c)+1)';
tout_c      = tout_c(1:end-1);
rbase_out_c = exp(1i*tout_c);

tin_c      = linspace(0, 2*pi, N_c+1)'; tin_c = tin_c(1:end-1);
tin_f      = linspace(0, 2*pi, N_f+1)'; tin_f = tin_f(1:end-1);
rbase_in_c = opt.Rp_c * exp(1i*tin_c);
rbase_in_f = opt.Rp_f * exp(1i*tin_f);

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q, opt);

% Assemble global grids
n_out      = numel(rbase_out_c);
rvec_in_c  = zeros(P*N_c, 1);
rout       = zeros(P*n_out, 1);
for k = 1:P
    rvec_in_c((k-1)*N_c+1:k*N_c)       = q(k) + rbase_in_c;
    rout((k-1)*n_out+1:k*n_out)         = q(k) + rbase_out_c;
end

svd_opts = struct(...
    'column_weight', logical(getOptField(opt,'column_weight',false)),...
    'left_weight',   logical(getOptField(opt,'left_weight',  false)));

[U, Y] = getSelfPseudo(1, rbase_in_c, rbase_out_c, [], [],...
    [0, n_out], 0, svd_opts);

opt_pair               = opt;
opt_pair.project_force = false;

[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU,pair_cache] =...
    getPairBasisStokes(q, rbase_in_c, rbase_in_f, rimage_vec, refine,...
    pairs, opt_pair, [], rbase_out_c, svd_opts);


geom            = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine     = refine;
geom.rimage_vec = rimage_vec;
geom.opt        = opt_pair;
geom.rvec_out   = rout;
geom.rcheck     = rout;
geom.q          = q;
geom.pairs      = pairs;
geom.rvec_in    = rvec_in_c;

basis            = struct();
basis.U          = U;
basis.Y          = Y;
basis.Upf        = UB_all;
basis.Ypf        = YB_all;
basis.DC_all     = UC_all;
basis.YC_all     = YC_all;
basis.Cmap       = Cmap;
basis.Cmap_FU    = Cmap_FU;
basis.pair_cache = pair_cache;

end

% -------------------------------------------------------------------------

function plotResBigSparseSpyFigure(big_sparse, stats, P, phi, case_info, plot_magn)
% Collect matrices and titles based on active mode
matrices = {};
titles   = {};

if case_info.sparse_map_coarse
    matrices{end+1} = big_sparse.M_source_corr;
    titles{end+1}   = 'Source correction map (direct coarse)';
end
matrices{end+1} = big_sparse.M_pair;
titles{end+1}   = 'Pair source map (used by velocity/FT correction)';

if strcmp(case_info.u_corr_mode,'combined')
    figure_name     = 'buildResPeanutBigSparseStokes: combined u';
    matrices{end+1} = big_sparse.M_u_corr;
    titles{end+1}   = 'Velocity correction map (combined)';
else
    figure_name     = 'buildResPeanutBigSparseStokes: factored u';
    matrices{end+1} = big_sparse.M_u_cross;
    titles{end+1}   = 'Velocity correction: cross part';
    matrices{end+1} = big_sparse.M_u_peanut;
    titles{end+1}   = 'Velocity correction: peanut part';
end

matrices{end+1} = big_sparse.M_ft_corr;
titles{end+1}   = 'Force/torque correction map';

[nrows, ncols] = chooseSpyTileShape(numel(matrices));
sup_lines = {...
    sprintf('P=%d, phi=%.2f, pairs=%d, N_c=%d', P, phi, stats.n_pairs, stats.N_c),...
    case_info.label};

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
    sup_lines{end+1} =...
        'WARNING: large magnitudes may indicate small N_c';
end
sgtitle(sup_lines);
end

function [nrows, ncols] = chooseSpyTileShape(nplots)
nrows = 1 + (nplots > 3);
ncols = ceil(nplots / nrows);
end
