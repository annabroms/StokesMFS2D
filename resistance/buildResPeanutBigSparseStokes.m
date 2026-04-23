function [big_sparse,stats] = buildResPeanutBigSparseStokes(geom,basis)
%BUILDRESPEANUTBIGSPARSESTOKES Build solve-grid resistance pair maps.
%
% Calling buildResPeanutBigSparseStokes with no inputs runs a small visual
% self-test and shows the sparse matrix structures with spy.

if nargin == 0
    test_buildResPeanutBigSparseStokes;
    big_sparse = [];
    stats = [];
    return
end

timer = tic;

% Validate the solve-grid contract before allocating sparse storage.
opt = geom.opt;
q = geom.q(:);
pairs = geom.pairs;
rcheck = geom.rcheck;
P = numel(q);
N_c = opt.N_c;
N_check = numel(rcheck)/P;
n_pairs = size(pairs,1);
n_coarse = P*N_c;
pair_rows = 4*N_c;
pair_total_rows = n_pairs*pair_rows;
n_source_cols = 2*n_coarse;
n_global_u_rows = 2*P*N_check;

validateBigSparseInputs(geom,basis,N_check);

% Estimate sparse storage and initialize bookkeeping.
ram_estimate = estimateResPeanutBigSparseRamStokes( ...
    P,N_c,N_check,n_pairs,opt);
stats = initBigSparseStats(n_pairs,N_c,N_check,ram_estimate);
stats.requested = true;
fprintf(['buildResPeanutBigSparseStokes: assembling solve-grid ', ...
    'sparse pair maps (P=%d, close pairs=%d, N_c=%d, ', ...
    'N_check=%d, velocity=%s, chunk pairs=%d).\n'], ...
    P,n_pairs,N_c,N_check,stats.velocity_correction,stats.chunk_pairs);

max_build_bytes = getOptField(opt,'res_big_sparse_max_build_bytes', ...
    getOptField(opt,'big_sparse_max_build_bytes',inf));
if stats.big_sparse_build_bytes > max_build_bytes
    error('buildResPeanutBigSparseStokes:BuildMemoryLimit', ...
        ['Estimated sparse-entry build memory %.3g bytes exceeds ', ...
         'the configured limit %.3g bytes.'], ...
        stats.big_sparse_build_bytes,max_build_bytes);
end

if n_pairs == 0
    big_sparse = emptyBigSparse(pair_total_rows,n_source_cols, ...
        n_global_u_rows,ram_estimate);
    stats.active = true;
    stats.reason = 'no_close_pairs';
    stats.build_time = toc(timer);
    big_sparse.stats = stats;
    return
end

% Assemble one dense pair contribution at a time into global sparse maps.
rout_base_c = rcheck(1:N_check) - q(1);
entries = initSparseBuilders(P,N_c,N_check,n_pairs,opt);
for row = 1:n_pairs
    i = pairs(row,1);
    j = pairs(row,2);
    C_pair = basis.Cmap{i,j};
    Cmap_FU = basis.Cmap_FU{i,j};
    [Ucross,Ecolloc] = buildStokesCoarsePairDense( ...
        q,geom.rbase_in_c(:),rout_base_c,pairs,row);
    entries = appendPairBlocks(entries,pairs,row,N_c,N_check,P, ...
        C_pair,Cmap_FU,Ucross,Ecolloc);
end
entries = flushAllSparseBuilders(entries);

% Store the fixed sparse operators used by transform and postprocessing.
big_sparse = struct();
big_sparse.matrix_plan = ram_estimate.matrix_plan;
big_sparse.ram_estimate = ram_estimate;
big_sparse.N_c = N_c;
big_sparse.pairs = pairs;
big_sparse.source_scatter_rows = buildSourceScatterRows(pairs,N_c,P);
big_sparse.velocity_scatter_rows = buildVelocityScatterRows( ...
    pairs,N_check,P);
big_sparse.ft_scatter_rows = buildForceTorqueScatterRows(pairs,P);
big_sparse.M_pair = entries.pair.S;
if strcmp(ram_estimate.matrix_plan.u_corr_mode,'combined')
    big_sparse.M_u_corr = entries.u.S;
else
    big_sparse.M_u_cross = entries.u_cross.S;
    big_sparse.M_u_peanut = entries.u_peanut.S;
end
big_sparse.M_ft_corr = entries.ft.S;

stats.nnz_pair = nnz(big_sparse.M_pair);
if strcmp(ram_estimate.matrix_plan.u_corr_mode,'combined')
    stats.nnz_u = nnz(big_sparse.M_u_corr);
else
    stats.nnz_u_cross = nnz(big_sparse.M_u_cross);
    stats.nnz_u_peanut = nnz(big_sparse.M_u_peanut);
end
stats.nnz_ft = nnz(big_sparse.M_ft_corr);
stats.nnz_source_scatter = numel(big_sparse.source_scatter_rows);
stats.nnz_velocity_scatter = numel(big_sparse.velocity_scatter_rows);
stats.nnz_ft_scatter = numel(big_sparse.ft_scatter_rows);
stats.active = true;
stats.reason = '';
stats.build_time = toc(timer);
big_sparse.stats = stats;
end

function validateBigSparseInputs(geom,basis,N_check)
opt = geom.opt;

if ~logical(getOptField(opt,'cmap',false))
    error('buildResPeanutBigSparseStokes:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.cmap=1.');
end
if ~logical(getOptField(opt,'self_correct',false))
    error('buildResPeanutBigSparseStokes:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.self_correct=1.');
end
if getOptField(opt,'N_peanut',0) <= 0
    error('buildResPeanutBigSparseStokes:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.N_peanut > 0.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('buildResPeanutBigSparseStokes:UnsupportedGrid', ...
        'Resistance big sparse pair maps are only built on the solve grid.');
end
if N_check ~= round(N_check)
    error('buildResPeanutBigSparseStokes:BadGridSize', ...
        'rcheck length must be divisible by the number of particles.');
end
if size(geom.pairs,1) > 0 && ...
        (~isfield(basis,'Cmap') || isempty(basis.Cmap) || ...
         ~isfield(basis,'Cmap_FU') || isempty(basis.Cmap_FU))
    error('buildResPeanutBigSparseStokes:MissingPairMaps', ...
        'Resistance big sparse assembly requires actual per-pair Cmap and Cmap_FU data.');
end
if size(geom.pairs,1) > 0 && isfield(basis,'pair_cache') && ...
        isfield(basis.pair_cache,'enabled') && basis.pair_cache.enabled
    error('buildResPeanutBigSparseStokes:PairCacheUnsupported', ...
        ['opt.use_big_sparse=1 requires actual per-pair maps. ', ...
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
stats.reason = 'not_requested';
stats.matrix_plan = plan;
stats.chunk_pairs = plan.chunk_pairs;
stats.n_pairs = n_pairs;
stats.N_c = N_c;
stats.N_check = N_check;
stats.used_pair_cache = false;
stats.source_correction = plan.source_correction;
stats.velocity_correction = plan.velocity_correction;
stats.u_corr_mode = plan.u_corr_mode;
stats.combined_u_corr = plan.combined_u_corr;
stats.force_torque_correction = plan.force_torque_correction;
stats.local_pair_entries = counts.pair/max(1,n_pairs);
stats.local_u_entries = counts.u/max(1,n_pairs);
stats.local_u_cross_entries = counts.u_cross/max(1,n_pairs);
stats.local_u_peanut_entries = counts.u_peanut/max(1,n_pairs);
stats.local_ft_entries = counts.ft/max(1,n_pairs);
stats.nnz_pair = 0;
stats.nnz_u = 0;
stats.nnz_u_cross = 0;
stats.nnz_u_peanut = 0;
stats.nnz_ft = 0;
stats.nnz_source_scatter = 0;
stats.nnz_velocity_scatter = 0;
stats.nnz_ft_scatter = 0;
stats.big_sparse_matrix_bytes = ram_estimate.big_sparse_matrix_bytes;
stats.big_sparse_auxiliary_bytes = ram_estimate.big_sparse_auxiliary_bytes;
stats.big_sparse_build_bytes = ram_estimate.big_sparse_build_bytes;
stats.big_sparse_peak_bytes = ram_estimate.big_sparse_peak_bytes;
stats.build_time = 0;
end

function big_sparse = emptyBigSparse(pair_total_rows,n_source_cols, ...
    n_global_u_rows,ram_estimate)
big_sparse = struct();
big_sparse.matrix_plan = ram_estimate.matrix_plan;
big_sparse.ram_estimate = ram_estimate;
big_sparse.N_c = [];
big_sparse.pairs = zeros(0,2);
big_sparse.source_scatter_rows = zeros(0,1);
big_sparse.velocity_scatter_rows = zeros(0,1);
big_sparse.ft_scatter_rows = zeros(0,1);
big_sparse.M_pair = sparse(pair_total_rows,n_source_cols);
switch ram_estimate.matrix_plan.u_corr_mode
    case 'combined'
        big_sparse.M_u_corr = sparse(0,n_source_cols + pair_total_rows);
    otherwise
        big_sparse.M_u_cross = sparse(0,n_source_cols);
        big_sparse.M_u_peanut = sparse(0,pair_total_rows);
end
big_sparse.M_ft_corr = sparse(0,n_source_cols);
big_sparse.n_global_u_rows = n_global_u_rows;
end

function entries = initSparseBuilders(P,N_c,N_check,n_pairs,opt)
chunk_pairs = max(1,round(getOptField(opt,'res_big_sparse_chunk_pairs',8)));
chunk_pairs = min(max(1,n_pairs),chunk_pairs);

n_coarse = P*N_c;
n_source_cols = 2*n_coarse;
pair_rows = 4*N_c;
pair_total_rows = n_pairs*pair_rows;
pair_u_rows = n_pairs*4*N_check;
pair_ft_rows = n_pairs*6;

local_pair = pair_rows*(4*N_c);
local_u = (4*N_check)*(4*N_c);
local_ft = 6*(4*N_c);
u_corr_mode = resolveResBigSparseUCorrMode(opt);

entries = struct();
entries.u_corr_mode = u_corr_mode;
entries.pair = initSparseBuilder( ...
    chunk_pairs*local_pair,pair_total_rows,n_source_cols);
switch u_corr_mode
    case 'combined'
        entries.u = initSparseBuilder( ...
            chunk_pairs*2*local_u,pair_u_rows,n_source_cols + pair_total_rows);
        entries.u_cross = initSparseBuilder(0,pair_u_rows,n_source_cols);
        entries.u_peanut = initSparseBuilder(0,pair_u_rows,pair_total_rows);
    otherwise
        entries.u = initSparseBuilder(0,pair_u_rows,n_source_cols);
        entries.u_cross = initSparseBuilder( ...
            chunk_pairs*local_u,pair_u_rows,n_source_cols);
        entries.u_peanut = initSparseBuilder( ...
            chunk_pairs*local_u,pair_u_rows,pair_total_rows);
end
entries.ft = initSparseBuilder( ...
    chunk_pairs*local_ft,pair_ft_rows,n_source_cols);
end

function entries = appendPairBlocks(entries,pairs,row,N_c,N_check,P, ...
    C_pair,Cmap_FU,Ucross,Ecolloc)
i = pairs(row,1);
j = pairs(row,2);
in_idx = pairCoarseInputIndices(i,j,N_c,P);
pair_idx = (row-1)*(4*N_c)+1:row*(4*N_c);
u_idx = pairVelocityBlockIndices(row,N_check);
ft_idx = pairForceTorqueBlockIndices(row);

entries.pair = appendSparseBuilderBlock(entries.pair, ...
    pair_idx,in_idx,C_pair);
switch entries.u_corr_mode
    case 'combined'
        n_source_cols = 2*P*N_c;
        entries.u = appendSparseBuilderBlock(entries.u,u_idx,in_idx, ...
            Ucross);
        entries.u = appendSparseBuilderBlock(entries.u,u_idx, ...
            n_source_cols + pair_idx,-Ecolloc);
    otherwise
        entries.u_cross = appendSparseBuilderBlock(entries.u_cross, ...
            u_idx,in_idx,Ucross);
        entries.u_peanut = appendSparseBuilderBlock(entries.u_peanut, ...
            u_idx,pair_idx,Ecolloc);
end
entries.ft = appendSparseBuilderBlock(entries.ft, ...
    ft_idx,in_idx,Cmap_FU);
end

function entries = flushAllSparseBuilders(entries)
entries.pair = flushSparseBuilderBlock(entries.pair);
entries.u = flushSparseBuilderBlock(entries.u);
entries.u_cross = flushSparseBuilderBlock(entries.u_cross);
entries.u_peanut = flushSparseBuilderBlock(entries.u_peanut);
entries.ft = flushSparseBuilderBlock(entries.ft);
end

function rows = buildSourceScatterRows(pairs,N_c,P)
n_pairs = size(pairs,1);
rows = zeros(n_pairs*4*N_c,1);
next = 1;
for row = 1:n_pairs
    idx = pairSourceOutputIndices(pairs(row,1),pairs(row,2),N_c,P);
    rows(next:next+numel(idx)-1) = idx;
    next = next + numel(idx);
end
end

function rows = buildVelocityScatterRows(pairs,N_check,P)
n_pairs = size(pairs,1);
rows = zeros(n_pairs*4*N_check,1);
next = 1;
for row = 1:n_pairs
    idx = pairVelocityOutputIndices(pairs(row,1),pairs(row,2), ...
        N_check,P);
    rows(next:next+numel(idx)-1) = idx;
    next = next + numel(idx);
end
end

function rows = buildForceTorqueScatterRows(pairs,P)
n_pairs = size(pairs,1);
rows = zeros(n_pairs*6,1);
next = 1;
for row = 1:n_pairs
    idx = pairForceTorqueOutputIndices(pairs(row,1),pairs(row,2),P);
    rows(next:next+numel(idx)-1) = idx;
    next = next + numel(idx);
end
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
idx = [idx_i, idx_j, n_coarse+idx_i, n_coarse+idx_j]';
end

function idx = pairVelocityOutputIndices(i,j,N_check,P)
pm_check = P*N_check;
idx_i = (i-1)*N_check+1:i*N_check;
idx_j = (j-1)*N_check+1:j*N_check;
idx = [idx_i, idx_j, pm_check+idx_i, pm_check+idx_j]';
end

function idx = pairForceTorqueOutputIndices(i,j,P)
idx_i = (i-1)*3+1:3*i;
idx_j = (j-1)*3+1:3*j;
idx = [idx_i idx_j]';
if max(idx) > 3*P
    error('buildResPeanutBigSparseStokes:BadForceTorqueIndex', ...
        'Internal force/torque output index out of range.');
end
end

function idx = pairVelocityBlockIndices(row,N_check)
idx = (row-1)*(4*N_check)+1:row*(4*N_check);
idx = idx(:);
end

function idx = pairForceTorqueBlockIndices(row)
idx = (row-1)*6+1:row*6;
idx = idx(:);
end

function test_buildResPeanutBigSparseStokes
fprintf('buildResPeanutBigSparseStokes self-test: spy plots\n');

P = 20;
delta = 0.1;
q = grow_cluster(P,delta,2);
opt = get2Dparams(P,18,24);
opt.delta_pair = 0.12;
opt.N_peanut = 60;
opt.get_bndry_field = 0;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.reuse_pair_basis_by_sep = false;
opt.use_fmm = false;

[geom,basis] = buildSelfTestData(q,opt);
[big_sparse,stats] = buildResPeanutBigSparseStokes(geom,basis);

fprintf('  P=%d, close pairs=%d\n',P,stats.n_pairs);
fprintf(['  M_pair nnz=%d, M_u_cross nnz=%d, ', ...
    'M_u_peanut nnz=%d, M_ft_corr nnz=%d\n'], ...
    stats.nnz_pair,stats.nnz_u_cross,stats.nnz_u_peanut, ...
    stats.nnz_ft);

figure('Name','buildResPeanutBigSparseStokes spy self-test','Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
spy(big_sparse.M_pair);
title('M_pair','Interpreter','none');

nexttile;
spy(big_sparse.M_u_cross);
title('M_u_cross','Interpreter','none');

nexttile;
spy(big_sparse.M_u_peanut);
title('M_u_peanut','Interpreter','none');

nexttile;
spy(big_sparse.M_ft_corr);
title('M_ft_corr','Interpreter','none');
end

function [geom,basis] = buildSelfTestData(q,opt)
q = q(:);
P = numel(q);
N_c = opt.N_c;
N_f = opt.N_f;
a_c = opt.a_c;

tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';
rbase_out_c = cos(tout_c)+1i*sin(tout_c);

tin_c = linspace(0,2*pi,N_c+1)';
tin_c = tin_c(1:end-1);
rbase_in_c = opt.Rp_c*cos(tin_c) + 1i*opt.Rp_c*sin(tin_c);

tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
rbase_in_f = opt.Rp_f*cos(tin_f) + 1i*opt.Rp_f*sin(tin_f);

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

rvec_in_c = zeros(P*N_c,1);
rout = zeros(P*numel(rbase_out_c),1);
for k = 1:P
    rvec_in_c((k-1)*N_c+1:k*N_c) = q(k) + rbase_in_c;
    out_idx = (k-1)*numel(rbase_out_c)+1:k*numel(rbase_out_c);
    rout(out_idx) = q(k) + rbase_out_c;
end

svd_opts = struct( ...
    'column_weight',logical(getOptField(opt,'column_weight',false)), ...
    'left_weight',logical(getOptField(opt,'left_weight',false)));
[U,Y] = getSelfPseudo(1,rbase_in_c,rbase_out_c,[],[], ...
    [0,numel(rbase_out_c)],0,svd_opts);

opt_pair = opt;
opt_pair.project_force = false;
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_FU,pair_cache] = ...
    getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine, ...
    pairs,opt_pair,[],rbase_out_c,svd_opts);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.rimage_vec = rimage_vec;
geom.opt = opt_pair;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.rvec_in = rvec_in_c;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Upf = UB_all;
basis.Ypf = YB_all;
basis.DC_all = UC_all;
basis.YC_all = YC_all;
basis.Cmap = Cmap;
basis.Cmap_FU = Cmap_FU;
basis.pair_cache = pair_cache;
end
