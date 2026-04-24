function [big_sparse,stats] = buildMobPeanutBigSparseStokes(geom,basis)
%BUILDMOBPEANUTBIGSPARSESTOKES Build solve-grid sparse pair-correction maps.
%
% Calling buildMobPeanutBigSparseStokes with no inputs runs a visual
% self-test on a small random_discs_mc packing at phi=0.65 and shows the
% sparse matrix structures with spy.
%
% This builder is deliberately actual-pair only. Rotations are not part of
% this path: solve_mob_peanut_enhanced disables
% opt.reuse_pair_basis_by_sep when opt.use_big_sparse is true, so every
% close pair contributes the dense maps for its actual geometry.
%
% opt.sparse_map_coarse and opt.big_sparse_direct_u_corr select the matrix
% set built here:
%   sparse_map_coarse = 0:
%       M_pair_nonp, P_pair, and source_scatter_rows keep the pair based
%       row structure
%       coarse correction
%           pair_nonp = M_pair_nonp*lambda_self
%           pair_proj = P_pair*pair_nonp
%           corr = scatter([pair_proj; pair_nonp]).
%   sparse_map_coarse = 1:
%       M_source_corr maps projected one-body coarse sources directly to
%       projected coarse-source increments, so the solve matvec applies
%           corr = M_source_corr*lambda_self.
%   big_sparse_direct_u_corr = 1:
%       M_u_corr maps projected one-body coarse sources directly to u_corr.
%   big_sparse_direct_u_corr = 0:
%       M_u_cross and M_u_peanut keep the velocity correction factored.
%       When sparse_map_coarse = 1 this also builds M_pair_proj so the
%       factored u_corr path can reuse projected pair data without a
%       runtime dense projector.
%   M_rbm_corr is always built because the solve-grid mobility path still
%       uses it when postprocessing rigid-body motion.
%
% No canonical-pair rotations are built in this path.

if nargin == 0
    test_buildMobPeanutBigSparseStokes;
    big_sparse = [];
    stats = [];
    return
end

timer = tic;

opt = geom.opt;
q = geom.q(:);
pairs = geom.pairs;
rcheck = geom.rcheck;
N_c = opt.N_c;
P = numel(q);
n_pairs = size(pairs,1);
n_coarse = P*N_c;
N_check = numel(rcheck)/P;
build_mode = resolveMobBigSparseBuildMode(opt);

validateBigSparseInputs(geom,basis,N_check,build_mode);

n_source_cols = 2*n_coarse;
n_u_rows = 2*P*N_check;
n_u_cols = 2*n_coarse;
pair_block_rows = 4*N_c;
pair_total_rows = n_pairs*pair_block_rows;
ram_estimate = estimateMobPeanutBigSparseRamStokes(P,N_c,N_check,n_pairs,opt);
plan = ram_estimate.matrix_plan;
direct_u_corr = plan.direct_u_corr;
sparse_map_coarse = plan.sparse_map_coarse;
need_pair_projection_matrix = n_pairs > 0;
if need_pair_projection_matrix
    P_pair = resolvePairProjectionMatrix(geom,basis);
else
    P_pair = [];
end

stats = initBigSparseStats(n_pairs,N_c,N_check,ram_estimate);
stats.requested = true;
fprintf(['buildMobPeanutBigSparseStokes: assembling solve-grid ', ...
    'sparse pair maps (mode=%s, P=%d, close pairs=%d, N_c=%d, ', ...
    'N_check=%d, chunk pairs=%d).\n'], ...
    build_mode,P,n_pairs,N_c,N_check,stats.chunk_pairs);

max_build_bytes = getOptField(opt,'big_sparse_max_build_bytes',inf);
if stats.big_sparse_build_bytes > max_build_bytes
    error('buildMobPeanutBigSparseStokes:BuildMemoryLimit', ...
        ['Estimated sparse-entry build memory %.3g bytes exceeds ', ...
         'opt.big_sparse_max_build_bytes %.3g bytes.'], ...
        stats.big_sparse_build_bytes,max_build_bytes);
end

if n_pairs == 0
    big_sparse = emptyBigSparse(n_source_cols,n_u_rows,n_u_cols, ...
        pair_total_rows,3*P,plan);
    big_sparse.ram_estimate = ram_estimate;
    stats.active = true;
    stats.reason = 'no_close_pairs';
    stats.build_time = toc(timer);
    big_sparse.stats = stats;
    return
end

rout_base_c = rcheck(1:N_check) - q(1);
entries = initSparseBuilders(P,N_c,N_check,n_pairs,opt,plan);

switch build_mode
    case 'precomputed'
        entries = appendPrecomputedBlocks(entries,geom,basis, ...
            rout_base_c,P_pair);
    case 'streaming'
        entries = appendStreamingBlocks(entries,geom,basis, ...
            rout_base_c,P_pair);
    otherwise
        error('buildMobPeanutBigSparseStokes:BadMode', ...
            'Unsupported build mode "%s".',build_mode);
end

entries = flushAllSparseBuilders(entries);

big_sparse = struct();
big_sparse.matrix_plan = plan;
big_sparse.ram_estimate = ram_estimate;
if sparse_map_coarse
    big_sparse.M_source_corr = entries.source_corr.S;
else
    big_sparse.P_pair = P_pair;
    big_sparse.source_scatter_rows = buildSourceScatterRows(pairs,N_c,P);
    big_sparse.M_pair_nonp = entries.pair_nonp.S;
end
big_sparse.M_rbm_corr = entries.rbm.S;
if direct_u_corr
    big_sparse.M_u_corr = entries.u.S;
else
    big_sparse.M_u_cross = entries.u_cross.S;
    big_sparse.M_u_peanut = entries.u_peanut.S;
    if sparse_map_coarse
        big_sparse.M_pair_proj = entries.pair_proj.S;
    end
end

stats.nnz_source_corr = 0;
stats.nnz_pair_nonp = 0;
stats.nnz_pair_proj = 0;
if direct_u_corr
    stats.nnz_u = nnz(big_sparse.M_u_corr);
else
    stats.nnz_u_cross = nnz(big_sparse.M_u_cross);
    stats.nnz_u_peanut = nnz(big_sparse.M_u_peanut);
    if sparse_map_coarse
        stats.nnz_pair_proj = nnz(big_sparse.M_pair_proj);
    end
end
if sparse_map_coarse
    stats.nnz_source_corr = nnz(big_sparse.M_source_corr);
else
    stats.nnz_pair_nonp = nnz(big_sparse.M_pair_nonp);
    stats.nnz_source_scatter = numel(big_sparse.source_scatter_rows);
end
stats.nnz_rbm = nnz(big_sparse.M_rbm_corr);
stats.active = true;
stats.reason = '';
stats.build_time = toc(timer);
big_sparse.stats = stats;
end

function build_mode = resolveMobBigSparseBuildMode(opt)
build_mode = lower(char(getOptField(opt,'mob_big_sparse_build_mode', ...
    'precomputed')));
switch build_mode
    case {'precomputed','streaming'}
        % supported
    otherwise
        error('buildMobPeanutBigSparseStokes:BadMode', ...
            ['opt.mob_big_sparse_build_mode must be ''precomputed'' ', ...
             'or ''streaming''.']);
end
end

function validateBigSparseInputs(geom,basis,N_check,build_mode)
opt = geom.opt;

if ~logical(getOptField(opt,'cmap',false))
    error('buildMobPeanutBigSparseStokes:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.cmap=1.');
end
if ~logical(getOptField(opt,'self_correct',false))
    error('buildMobPeanutBigSparseStokes:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.self_correct=1.');
end
if ~logical(getOptField(opt,'use_dense',false))
    error('buildMobPeanutBigSparseStokes:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.use_dense=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('buildMobPeanutBigSparseStokes:UnsupportedGrid', ...
        'opt.use_big_sparse=1 is only supported on the solve grid.');
end
if N_check ~= round(N_check)
    error('buildMobPeanutBigSparseStokes:BadGridSize', ...
        'rcheck length must be divisible by the number of particles.');
end
if strcmp(build_mode,'streaming') && ...
        logical(getOptField(opt,'get_bndry_field',false))
    error('buildMobPeanutBigSparseStokes:StreamingBoundaryUnsupported', ...
        ['opt.mob_big_sparse_build_mode=''streaming'' requires ', ...
         'opt.get_bndry_field=0.']);
end
if strcmp(build_mode,'streaming') && ...
        logical(getOptField(opt,'reuse_pair_basis_by_sep',false))
    error('buildMobPeanutBigSparseStokes:StreamingPairReuseUnsupported', ...
        ['opt.mob_big_sparse_build_mode=''streaming'' requires ', ...
         'opt.reuse_pair_basis_by_sep=0.']);
end
if strcmp(build_mode,'precomputed') && size(geom.pairs,1) > 0 && ...
        (~isfield(basis,'Cmap') || isempty(basis.Cmap) || ...
         ~isfield(basis,'Cmap_FU') || isempty(basis.Cmap_FU))
    error('buildMobPeanutBigSparseStokes:MissingCmap', ...
        'Precomputed big sparse requires actual per-pair Cmap and Cmap_FU data.');
end
if strcmp(build_mode,'streaming') && size(geom.pairs,1) > 0
    required_geom = {'rbase_in_f','rimage_vec','refine'};
    for k = 1:numel(required_geom)
        if ~isfield(geom,required_geom{k}) || isempty(geom.(required_geom{k}))
            error('buildMobPeanutBigSparseStokes:MissingStreamingGeom', ...
                'Streaming big sparse requires geom.%s.',required_geom{k});
        end
    end
    if ~isfield(basis,'Lc') || isempty(basis.Lc)
        error('buildMobPeanutBigSparseStokes:MissingLc', ...
            'Streaming big sparse requires basis.Lc.');
    end
end
if size(geom.pairs,1) > 0 && isfield(basis,'pair_cache') && ...
        isfield(basis.pair_cache,'enabled') && basis.pair_cache.enabled
    error('buildMobPeanutBigSparseStokes:PairCacheUnsupported', ...
        ['opt.use_big_sparse=1 currently requires actual per-pair maps. ', ...
         'Set opt.reuse_pair_basis_by_sep=0 before building the pair basis.']);
end
end

function stats = initBigSparseStats(n_pairs,N_c,N_check,ram_estimate)
counts = ram_estimate.counts;
plan = ram_estimate.matrix_plan;
stats = struct();
stats.requested = false;
stats.active = false;
stats.backend = 'global_block_sparse';
stats.matrix_plan = plan;
stats.build_mode = plan.build_mode;
stats.chunk_pairs = plan.chunk_pairs;
stats.reason = 'not_requested';
stats.n_pairs = n_pairs;
stats.N_c = N_c;
stats.N_check = N_check;
stats.used_pair_cache = false;
stats.source_correction = plan.source_correction;
stats.velocity_correction = plan.velocity_correction;
stats.direct_u_corr = plan.direct_u_corr;
stats.sparse_map_coarse = plan.sparse_map_coarse;
stats.projector_mode = plan.projector_mode;
stats.local_source_corr_entries = counts.source_corr/max(1,n_pairs);
stats.local_pair_nonp_entries = counts.pair_nonp/max(1,n_pairs);
stats.local_pair_proj_entries = counts.pair_proj/max(1,n_pairs);
stats.local_u_entries = counts.u/max(1,n_pairs);
stats.local_u_cross_entries = counts.u_cross/max(1,n_pairs);
stats.local_u_peanut_entries = counts.u_peanut/max(1,n_pairs);
stats.local_rbm_entries = counts.rbm/max(1,n_pairs);
stats.nnz_source_corr = 0;
stats.nnz_u = 0;
stats.nnz_u_cross = 0;
stats.nnz_u_peanut = 0;
stats.nnz_pair_nonp = 0;
stats.nnz_pair_proj = 0;
stats.nnz_rbm = 0;
stats.nnz_source_scatter = 0;
stats.rotations_used = false;
stats.big_sparse_matrix_bytes = ram_estimate.big_sparse_matrix_bytes;
stats.big_sparse_auxiliary_bytes = ram_estimate.big_sparse_auxiliary_bytes;
stats.big_sparse_build_bytes = ram_estimate.big_sparse_build_bytes;
stats.big_sparse_peak_bytes = ram_estimate.big_sparse_peak_bytes;
stats.retained_pair_basis_bytes = ram_estimate.retained_pair_basis_bytes;
stats.solver_precompute_peak_bytes = ...
    ram_estimate.solver_precompute_peak_bytes;
stats.build_time = 0;
end

function big_sparse = emptyBigSparse(n_source_cols,n_u_rows,n_u_cols, ...
    pair_total_rows,n_rbm_rows,plan)
big_sparse = struct();
big_sparse.matrix_plan = plan;
if plan.sparse_map_coarse
    big_sparse.M_source_corr = sparse(n_source_cols,n_source_cols);
else
    big_sparse.P_pair = sparse(0,0);
    big_sparse.source_scatter_rows = zeros(0,1);
    big_sparse.M_pair_nonp = sparse(pair_total_rows,n_source_cols);
end
big_sparse.M_rbm_corr = sparse(n_rbm_rows,n_source_cols);
if plan.direct_u_corr
    big_sparse.M_u_corr = sparse(n_u_rows,n_u_cols);
else
    big_sparse.M_u_cross = sparse(n_u_rows,n_u_cols);
    big_sparse.M_u_peanut = sparse(n_u_rows,pair_total_rows);
    if plan.sparse_map_coarse
        big_sparse.M_pair_proj = sparse(pair_total_rows,n_source_cols);
    end
end
end

function P_pair = resolvePairProjectionMatrix(geom,basis)
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
    error('buildMobPeanutBigSparseStokes:MissingDensePairProjector', ...
        ['opt.use_matrix_free_Lc_pair=0 requires basis.Lc_pair ', ...
         'or basis.Lc.']);
end

expected = 4*geom.opt.N_c;
if ~isequal(size(P_pair),[expected expected])
    error('buildMobPeanutBigSparseStokes:BadDensePairProjector', ...
        'basis.Lc_pair must have size %d-by-%d.',expected,expected);
end
end

function entries = initSparseBuilders(P,N_c,N_check,n_pairs,opt,plan)
chunk_pairs = max(1,round(getOptField(opt,'mob_big_sparse_chunk_pairs',8)));
chunk_pairs = min(max(1,n_pairs),chunk_pairs);
n_coarse = P*N_c;
n_source_cols = 2*n_coarse;
n_u_rows = 2*P*N_check;
n_u_cols = 2*n_coarse;
pair_rows = 4*N_c;
pair_total_rows = n_pairs*pair_rows;
n_rbm_rows = 3*P;
local_pair_nonp = pair_rows*(4*N_c);
local_u = (4*N_check)*(4*N_c);
local_rbm = 6*(4*N_c);

entries = struct();
if plan.sparse_map_coarse
    entries.source_corr = initSparseBuilder( ...
        chunk_pairs*local_pair_nonp,n_source_cols,n_source_cols);
    entries.pair_nonp = initSparseBuilder(0,pair_total_rows,n_source_cols);
else
    entries.source_corr = initSparseBuilder(0,n_source_cols,n_source_cols);
    entries.pair_nonp = initSparseBuilder( ...
        chunk_pairs*local_pair_nonp,pair_total_rows,n_source_cols);
end
entries.rbm = initSparseBuilder( ...
    chunk_pairs*local_rbm,n_rbm_rows,n_source_cols);
if plan.direct_u_corr
    entries.u = initSparseBuilder(chunk_pairs*local_u,n_u_rows,n_u_cols);
    entries.u_cross = initSparseBuilder(0,n_u_rows,n_u_cols);
    entries.u_peanut = initSparseBuilder(0,n_u_rows,pair_total_rows);
    entries.pair_proj = initSparseBuilder(0,pair_total_rows,n_source_cols);
else
    entries.u = initSparseBuilder(0,n_u_rows,n_u_cols);
    entries.u_cross = initSparseBuilder(chunk_pairs*local_u,n_u_rows,n_u_cols);
    entries.u_peanut = initSparseBuilder(chunk_pairs*local_u,n_u_rows, ...
        pair_total_rows);
    if plan.sparse_map_coarse
        entries.pair_proj = initSparseBuilder( ...
            chunk_pairs*local_pair_nonp,pair_total_rows,n_source_cols);
    else
        entries.pair_proj = initSparseBuilder(0,pair_total_rows,n_source_cols);
    end
end
end

function entries = appendPrecomputedBlocks(entries,geom,basis, ...
    rout_base_c,P_pair)
pairs = geom.pairs;
q = geom.q(:);
N_c = geom.opt.N_c;
N_check = numel(geom.rcheck)/numel(q);
P = numel(q);
n_pairs = size(pairs,1);

for row = 1:n_pairs
    i = pairs(row,1);
    j = pairs(row,2);
    C_nonp = basis.Cmap{i,j};
    Cmap_FU = basis.Cmap_FU{i,j};
    [Ucross,Ecolloc] = buildStokesCoarsePairDense( ...
        q,geom.rbase_in_c(:),rout_base_c,pairs,row);
    entries = appendPairBlocks(entries,pairs,row,N_c,N_check,P, ...
        C_nonp,Cmap_FU,Ucross,Ecolloc,P_pair,geom.opt);
end
end

function entries = appendStreamingBlocks(entries,geom,basis, ...
    rout_base_c,P_pair)
pairs = geom.pairs;
q = geom.q(:);
opt = geom.opt;
N_c = opt.N_c;
N_check = numel(geom.rcheck)/numel(q);
P = numel(q);
n_pairs = size(pairs,1);
show_counter = logical(getOptField(opt,'show_counter',false));
svd_opts = struct( ...
    'column_weight',logical(getOptField(opt,'column_weight',false)), ...
    'left_weight',logical(getOptField(opt,'left_weight',false)));

pair_opt = opt;
pair_opt.project_force = true;
pair_opt.project = true;
pair_opt.pair_basis_debug = false;

for row = 1:n_pairs
    pair = buildStokesMobilityPairData(q,geom.rbase_in_c(:), ...
        geom.rbase_in_f(:),geom.rimage_vec,geom.refine,pairs,pair_opt, ...
        basis.Lc,row,false,svd_opts,false,'maps_only');
    [Ucross,Ecolloc] = buildStokesCoarsePairDense( ...
        q,geom.rbase_in_c(:),rout_base_c,pairs,row);
    entries = appendPairBlocks(entries,pairs,row,N_c,N_check,P, ...
        pair.Cmap,pair.Cmap_FU,Ucross,Ecolloc,P_pair,opt);
    if show_counter
        fprintf('buildMobPeanutBigSparseStokes: streamed pair %d/%d\n', ...
            row,n_pairs);
    end
end
end

function entries = appendPairBlocks(entries,pairs,row,N_c,N_check,P, ...
    C_nonp,Cmap_FU,Ucross,Ecolloc,P_pair,opt)
i = pairs(row,1);
j = pairs(row,2);
in_idx = pairCoarseInputIndices(i,j,N_c,P);
proj_out_idx = pairProjectedSourceOutputIndices(i,j,N_c,P);
u_out_idx = pairVelocityOutputIndices(i,j,N_check,P);
pair_idx = (row-1)*(4*N_c)+1:row*(4*N_c);
rbm_idx = pairRigidOutputIndices(i,j,P);
direct_u_corr = logical(getOptField(opt,'big_sparse_direct_u_corr',true));
sparse_map_coarse = logical(getOptField(opt,'sparse_map_coarse',false));

if sparse_map_coarse || direct_u_corr
    C_proj = P_pair*C_nonp;
else
    C_proj = [];
end
if sparse_map_coarse
    entries.source_corr = appendSparseBuilderBlock(entries.source_corr, ...
        proj_out_idx,in_idx,C_proj);
else
    entries.pair_nonp = appendSparseBuilderBlock(entries.pair_nonp, ...
        pair_idx,in_idx,C_nonp);
end
entries.rbm = appendSparseBuilderBlock(entries.rbm,rbm_idx,in_idx,-Cmap_FU);
if direct_u_corr
    entries.u = appendSparseBuilderBlock(entries.u,u_out_idx,in_idx, ...
        Ucross - Ecolloc*C_proj);
else
    entries.u_cross = appendSparseBuilderBlock(entries.u_cross,u_out_idx,in_idx, ...
        Ucross);
    entries.u_peanut = appendSparseBuilderBlock(entries.u_peanut,u_out_idx, ...
        pair_idx,Ecolloc);
    if sparse_map_coarse
        entries.pair_proj = appendSparseBuilderBlock(entries.pair_proj, ...
            pair_idx,in_idx,C_proj);
    end
end
end

function entries = flushAllSparseBuilders(entries)
entries.source_corr = flushSparseBuilderBlock(entries.source_corr);
entries.pair_nonp = flushSparseBuilderBlock(entries.pair_nonp);
entries.rbm = flushSparseBuilderBlock(entries.rbm);
entries.u = flushSparseBuilderBlock(entries.u);
entries.u_cross = flushSparseBuilderBlock(entries.u_cross);
entries.u_peanut = flushSparseBuilderBlock(entries.u_peanut);
entries.pair_proj = flushSparseBuilderBlock(entries.pair_proj);
end

function rows = buildSourceScatterRows(pairs,N_c,P)
n_pairs = size(pairs,1);
rows = zeros(n_pairs*8*N_c,1);
next = 1;
for row = 1:n_pairs
    ij = pairSourceOutputIndices(pairs(row,1),pairs(row,2),N_c,P);
    rows(next:next+numel(ij)-1) = ij;
    next = next + numel(ij);
end
end

function P_pair = buildPairProjectionMatrix(rbase)
K = getKmat2D(rbase,0);
L = K*((K'*K)\K');
P_pair = getILpair(L);
end

function idx = pairCoarseInputIndices(i,j,N_c,P)
n_coarse = P*N_c;
idx_i = (i-1)*N_c+1:i*N_c;
idx_j = (j-1)*N_c+1:j*N_c;
idx = [idx_i, idx_j, n_coarse+idx_i, n_coarse+idx_j]';
end

function idx = pairSourceOutputIndices(i,j,N_c,P)
n_coarse = P*N_c;
idx_i = (i-1)*N_c+1:i*N_c;
idx_j = (j-1)*N_c+1:j*N_c;
proj_idx = [idx_i, idx_j, n_coarse+idx_i, n_coarse+idx_j];
nonp_idx = [2*n_coarse+idx_i, 2*n_coarse+idx_j, ...
            3*n_coarse+idx_i, 3*n_coarse+idx_j];
idx = [proj_idx, nonp_idx]';
end

function idx = pairProjectedSourceOutputIndices(i,j,N_c,P)
n_coarse = P*N_c;
idx_i = (i-1)*N_c+1:i*N_c;
idx_j = (j-1)*N_c+1:j*N_c;
idx = [idx_i, idx_j, n_coarse+idx_i, n_coarse+idx_j]';
end

function idx = pairVelocityOutputIndices(i,j,N_check,P)
pm_check = P*N_check;
idx_i = (i-1)*N_check+1:i*N_check;
idx_j = (j-1)*N_check+1:j*N_check;
idx = [idx_i, idx_j, pm_check+idx_i, pm_check+idx_j]';
end

function idx = pairRigidOutputIndices(i,j,P)
idx_i = (i-1)*3+1:3*i;
idx_j = (j-1)*3+1:3*j;
idx = [idx_i idx_j]';
if max(idx) > 3*P
    error('buildMobPeanutBigSparseStokes:BadRigidIndex', ...
        'Internal rigid-motion output index out of range.');
end
end

function test_buildMobPeanutBigSparseStokes
fprintf('buildMobPeanutBigSparseStokes self-test: random_discs_mc spy plots\n');

P = 20;
phi = 0.65;
rad = 1;
geom_opt = struct('domain','boxed','phi',phi,'rad',rad, ...
    'min_gap',1e-3,'n_sweeps',200,'rng_seed',210421, ...
    'visualise',false);
[q,meta] = random_discs_mc(P,geom_opt);

N_c = 24;
N_f = 60;
opt = get2Dparams(P,N_c,N_f);
opt.N_peanut = 120;
opt.delta_pair = 0.2;
opt.get_bndry_field = 0;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.use_matrix_free_Lc_pair = 1;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.use_big_sparse = true;

[geom,basis] = buildSelfTestMobilityData(q,opt);
cases = struct( ...
    'sparse_map_coarse',{false,false,true,true}, ...
    'direct_u_corr',{true,false,true,false}, ...
    'label',{ ...
        'pair based coarse row structure, direct u', ...
        'pair based coarse row structure, factored u', ...
        'direct coarse, direct u: SMALLEST storage and FASTEST apply', ...
        'direct coarse, factored u'});

fprintf('  P=%d, target phi=%.3f, actual phi=%.3f\n',P,phi,meta.phi);
for case_id = 1:numel(cases)
    geom_case = geom;
    geom_case.opt.sparse_map_coarse = cases(case_id).sparse_map_coarse;
    geom_case.opt.big_sparse_direct_u_corr = cases(case_id).direct_u_corr;
    [big_sparse,stats] = buildMobPeanutBigSparseStokes(geom_case,basis);
    fprintf(['  case %d: sparse_map_coarse=%d direct_u_corr=%d ', ...
        'pairs=%d source_corr nnz=%d pair_nonp nnz=%d pair_proj nnz=%d ', ...
        'u_direct nnz=%d u_cross nnz=%d u_peanut nnz=%d rbm nnz=%d\n'], ...
        case_id,cases(case_id).sparse_map_coarse, ...
        cases(case_id).direct_u_corr,stats.n_pairs, ...
        stats.nnz_source_corr,stats.nnz_pair_nonp,stats.nnz_pair_proj, ...
        stats.nnz_u,stats.nnz_u_cross,stats.nnz_u_peanut,stats.nnz_rbm);
    plotMobilityBigSparseSpyFigure(big_sparse,stats,P,meta.phi,cases(case_id));
end

end

function [geom,basis] = buildSelfTestMobilityData(q,opt)
q = q(:);
P = numel(q);
rad = 1;
N_c = opt.N_c;
N_f = opt.N_f;
a_c = opt.a_c;
Rp_c = opt.Rp_c;
Rp_f = opt.Rp_f;

tout_c_all = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c_all(1:end-1)';
rbase_out_c = rad*(cos(tout_c)+1i*sin(tout_c));
tin_c = linspace(0,2*pi,N_c+1)';
tin_c = tin_c(1:end-1);
rbase_in_c = Rp_c*cos(tin_c) + 1i*Rp_c*sin(tin_c);
tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
rbase_in_f = Rp_f*cos(tin_f) + 1i*Rp_f*sin(tin_f);

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

rvec_in_c = zeros(P*N_c,1);
rout = zeros(P*numel(rbase_out_c),1);
for k = 1:P
    rvec_in_c((k-1)*N_c+1:k*N_c) = q(k) + rbase_in_c;
    out_idx = (k-1)*numel(rbase_out_c)+1:k*numel(rbase_out_c);
    rout(out_idx) = q(k) + rbase_out_c;
end

svd_opts = struct('column_weight',false,'left_weight',false);
[U,Y,Lc] = getSelfPseudoMobilityStokes( ...
    1,q,rbase_in_c,rbase_out_c,[],[0,ceil(a_c*N_c)],svd_opts);

opt_pair = opt;
opt_pair.project_force = true;
opt_pair.project = true;
opt_pair.pair_basis_debug = 0;
opt_pair.rad = ones(P,1);
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU,pair_cache] = ...
    getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine, ...
    pairs,opt_pair,Lc{1},rbase_out_c,svd_opts);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.rvec_in = rvec_in_c;
geom.rimage_vec = rimage_vec;
geom.opt = opt_pair;
geom.opt.get_bndry_field = 0;
geom.opt.parallel_solve = false;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.refine = refine;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Lc = Lc{1};
basis.Upf = UB_all;
basis.Ypf = YB_all;
basis.DC_all = UC_all;
basis.YC_all = YC_all;
basis.Cmap = Cmap;
basis.Cmap_FU = Cmap_FU;
basis.Lc_pair = getILpair(Lc{1});
basis.pair_cache = pair_cache;
end

function plotMobilityBigSparseSpyFigure(big_sparse,stats,P,phi,case_info)
source_is_direct = case_info.sparse_map_coarse;
if source_is_direct
    matrices = {big_sparse.M_source_corr};
    titles = {'Source correction map'};
    source_note = 'source: M_source_corr';
else
    matrices = {big_sparse.M_pair_nonp};
    titles = {'Source correction map'};
    source_note = 'source: M_pair_nonp with P_pair and source_scatter_rows';
end

if case_info.direct_u_corr
    figure_name = 'buildMobPeanutBigSparseStokes spy direct maps';
    matrices{end+1} = big_sparse.M_u_corr;
    titles{end+1} = 'Velocity correction map';
    velocity_note = 'velocity: M_u_corr';
else
    figure_name = 'buildMobPeanutBigSparseStokes spy factored maps';
    if source_is_direct
        matrices{end+1} = big_sparse.M_pair_proj;
        titles{end+1} = 'Projected pair map for factored velocity correction';
    end
    matrices{end+1} = big_sparse.M_u_cross;
    titles{end+1} = 'factored: Part 1 of velocity correction map';
    matrices{end+1} = big_sparse.M_u_peanut;
    titles{end+1} = 'factored: Part 2 of velocity correction map';
    velocity_note = 'velocity: M_u_cross and M_u_peanut';
end

[nrows,ncols] = chooseSpyTileShape(numel(matrices));
figure('Name',sprintf('%s (%s)',figure_name,case_info.label), ...
    'Color','w');
tiledlayout(nrows,ncols,'TileSpacing','compact','Padding','compact');
for k = 1:numel(matrices)
    nexttile;
    spy(matrices{k});
    title(titles{k},'Interpreter','none');
end

sgtitle({ ...
    sprintf('P=%d, phi=%.2f, pairs=%d',P,phi,stats.n_pairs), ...
    case_info.label, ...
    sprintf('not shown: RBM correction')});
end

function [nrows,ncols] = chooseSpyTileShape(nplots)
if nplots <= 3
    nrows = 1;
else
    nrows = 2;
end
ncols = ceil(nplots/nrows);
end
