function [big_sparse,stats] = buildMobPeanutBigSparseStokes(geom,basis)
%BUILDMOBPEANUTBIGSPARSESTOKES Build solve-grid sparse pair-correction maps.
%
% This v1 builder is deliberately actual-pair only. Rotations are not part
% of this path: solve_mob_peanut_enhanced disables
% opt.reuse_pair_basis_by_sep when opt.use_big_sparse is true, so every
% close pair contributes the dense maps for its actual geometry.
%
% Optional direct maps are:
%   M_source_corr: [lambda_self_x; lambda_self_y] ->
%       [projected dx; projected dy; unprojected dx; unprojected dy], only
%       when opt.big_sparse_direct_source_corr is true.
%   M_u_corr: [lambda_self_x; lambda_self_y] -> u_corr, only when
%       opt.big_sparse_direct_u_corr is true.
%
% The default sparse factored matvec uses the maps below. With
% opt.big_sparse_structured_apply, M_pair_projector is replaced by the dense
% local P_pair block and M_source_scatter is replaced by source_scatter_rows.
% With opt.big_sparse_combine_pair_u_corr, M_pair_nonp and M_u_corr are
% replaced by one stacked matrix M_pair_nonp_u_corr for timing comparison.
% The builder only creates the matrices needed by the selected modes:
%   M_pair_nonp, M_pair_projector, M_source_scatter, M_u_cross, M_u_peanut.

timer = tic;

opt = geom.opt;
q = geom.q(:);
pairs = geom.pairs;
rbase_in_c = geom.rbase_in_c(:);
rcheck = geom.rcheck;
N_c = opt.N_c;
P = numel(q);
n_pairs = size(pairs,1);
n_coarse = P*N_c;
N_check = numel(rcheck)/P;

validateBigSparseInputs(geom,basis,N_check);

n_source_rows = 4*n_coarse;
n_source_cols = 2*n_coarse;
n_u_rows = 2*P*N_check;
n_u_cols = 2*n_coarse;
pair_block_rows = 4*N_c;
pair_total_rows = n_pairs*pair_block_rows;
P_pair = buildPairProjectionMatrix(rbase_in_c);
ram_estimate = estimateMobPeanutBigSparseRamStokes(P,N_c,N_check,n_pairs,opt);
plan = ram_estimate.matrix_plan;
counts = ram_estimate.counts;

stats = initBigSparseStats(n_pairs,N_c,N_check,ram_estimate);
stats.requested = true;

max_build_bytes = getOptField(opt,'big_sparse_max_build_bytes',inf);
if stats.estimated_build_bytes > max_build_bytes
    error('buildMobPeanutBigSparseStokes:BuildMemoryLimit', ...
        ['Estimated sparse-entry build memory %.3g bytes exceeds ', ...
         'opt.big_sparse_max_build_bytes %.3g bytes.'], ...
        stats.estimated_build_bytes,max_build_bytes);
end

if n_pairs == 0
    big_sparse = emptyBigSparse(n_source_rows,n_source_cols,n_u_rows, ...
        n_u_cols,pair_total_rows,plan);
    big_sparse.ram_estimate = ram_estimate;
    stats.active = true;
    stats.reason = 'no_close_pairs';
    stats.build_time = toc(timer);
    big_sparse.stats = stats;
    return
end

rout_base_c = rcheck(1:N_check) - q(1);
entries = preallocateSparseEntries(counts);

for row = 1:n_pairs
    i = pairs(row,1);
    j = pairs(row,2);

    in_idx = pairCoarseInputIndices(i,j,N_c,P);
    source_out_idx = pairSourceOutputIndices(i,j,N_c,P);
    u_out_idx = pairVelocityOutputIndices(i,j,N_check,P);
    pair_idx = (row-1)*pair_block_rows+1:row*pair_block_rows;

    C_nonp = basis.Cmap{i,j};
    if plan.direct_source_corr || plan.direct_u_corr
        C_proj = P_pair*C_nonp;
    end
    [Ucross,Ecolloc] = buildActualCoarsePairDense( ...
        q,rbase_in_c,rout_base_c,pairs,row);

    % opt.big_sparse_direct_source_corr uses this direct source map.
    % It replaces M_pair_nonp, M_pair_projector, and M_source_scatter for
    % the source update, but the pair maps may still be needed by factored
    % u_corr unless opt.big_sparse_direct_u_corr is also true.
    if plan.direct_source_corr
        entries = appendDenseBlock(entries,'source',source_out_idx,in_idx, ...
            [C_proj; C_nonp]);
    end

    % opt.big_sparse_direct_u_corr uses this direct velocity map.
    % This is mainly diagnostic; the default factored u path avoids forming
    % Ucross - Ecolloc*C_proj explicitly.
    if plan.direct_u_corr
        Ucorr = Ucross - Ecolloc*C_proj;
        if plan.combined_pair_u_corr
            entries = appendDenseBlock(entries,'pair_nonp_u_corr', ...
                pair_total_rows+u_out_idx,in_idx,Ucorr);
        else
            entries = appendDenseBlock(entries,'u',u_out_idx,in_idx,Ucorr);
        end
    end

    % Factored source and factored u paths share these pair source maps.
    if plan.combined_pair_u_corr
        entries = appendDenseBlock(entries,'pair_nonp_u_corr',pair_idx, ...
            in_idx,C_nonp);
    elseif plan.need_pair_nonp_matrix
        entries = appendDenseBlock(entries,'pair_nonp',pair_idx,in_idx, ...
            C_nonp);
    end
    if plan.need_pair_projector
        entries = appendDenseBlock(entries,'pair_projector',pair_idx, ...
            pair_idx,P_pair);
    end

    % opt.big_sparse_direct_source_corr=false uses this scatter to
    % convert [pair_proj; pair_nonp] into global source increments.
    if plan.need_source_scatter
        entries = appendSourceScatter(entries,source_out_idx,pair_idx, ...
            pair_total_rows);
    end

    % opt.big_sparse_direct_u_corr=false uses these maps for the
    % stable factored velocity correction.
    if plan.need_u_cross
        entries = appendDenseBlock(entries,'u_cross',u_out_idx,in_idx, ...
            Ucross);
    end
    if plan.need_u_peanut
        entries = appendDenseBlock(entries,'u_peanut',u_out_idx,pair_idx, ...
            Ecolloc);
    end
end

big_sparse = struct();
big_sparse.matrix_plan = plan;
big_sparse.ram_estimate = ram_estimate;
if plan.need_dense_pair_projector
    big_sparse.P_pair = P_pair;
end
if plan.need_source_scatter_indices
    big_sparse.source_scatter_rows = buildSourceScatterRows(pairs,N_c,P);
end
if plan.direct_source_corr
    big_sparse.M_source_corr = sparseFromEntries(entries.source, ...
        n_source_rows,n_source_cols);
end
if plan.direct_u_corr
    if plan.combined_pair_u_corr
        big_sparse.M_pair_nonp_u_corr = sparseFromEntries( ...
            entries.pair_nonp_u_corr,pair_total_rows+n_u_rows,n_u_cols);
        big_sparse.pair_nonp_rows = pair_total_rows;
    else
        big_sparse.M_u_corr = sparseFromEntries(entries.u,n_u_rows,n_u_cols);
    end
end
if plan.need_pair_nonp_matrix
    big_sparse.M_pair_nonp = sparseFromEntries(entries.pair_nonp, ...
        pair_total_rows,n_source_cols);
end
if plan.need_pair_projector
    big_sparse.M_pair_projector = sparseFromEntries(entries.pair_projector, ...
        pair_total_rows,pair_total_rows);
end
if plan.need_source_scatter
    big_sparse.M_source_scatter = sparseFromEntries(entries.source_scatter, ...
        n_source_rows,2*pair_total_rows);
end
if plan.need_u_cross
    big_sparse.M_u_cross = sparseFromEntries(entries.u_cross,n_u_rows, ...
        n_u_cols);
end
if plan.need_u_peanut
    big_sparse.M_u_peanut = sparseFromEntries(entries.u_peanut,n_u_rows, ...
        pair_total_rows);
end

if plan.direct_source_corr
    stats.nnz_source = nnz(big_sparse.M_source_corr);
end
if plan.direct_u_corr
    if plan.combined_pair_u_corr
        stats.nnz_pair_nonp_u_corr = nnz(big_sparse.M_pair_nonp_u_corr);
    else
        stats.nnz_u = nnz(big_sparse.M_u_corr);
    end
end
if plan.need_pair_nonp_matrix
    stats.nnz_pair_nonp = nnz(big_sparse.M_pair_nonp);
end
if plan.need_pair_projector
    stats.nnz_pair_projector = nnz(big_sparse.M_pair_projector);
end
if plan.need_source_scatter
    stats.nnz_source_scatter = nnz(big_sparse.M_source_scatter);
end
if plan.need_source_scatter_indices
    stats.nnz_source_scatter = numel(big_sparse.source_scatter_rows);
end
if plan.need_u_cross
    stats.nnz_u_cross = nnz(big_sparse.M_u_cross);
end
if plan.need_u_peanut
    stats.nnz_u_peanut = nnz(big_sparse.M_u_peanut);
end
stats.active = true;
stats.reason = '';
stats.build_time = toc(timer);
big_sparse.stats = stats;
end

function validateBigSparseInputs(geom,basis,N_check)
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
if ~logical(getOptField(opt,'use_matrix_free_Lc_pair',true))
    error('buildMobPeanutBigSparseStokes:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.use_matrix_free_Lc_pair=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('buildMobPeanutBigSparseStokes:UnsupportedGrid', ...
        'opt.use_big_sparse=1 is only supported on the solve grid.');
end
if N_check ~= round(N_check)
    error('buildMobPeanutBigSparseStokes:BadGridSize', ...
        'rcheck length must be divisible by the number of particles.');
end
if size(geom.pairs,1) > 0 && ...
        (~isfield(basis,'Cmap') || isempty(basis.Cmap))
    error('buildMobPeanutBigSparseStokes:MissingCmap', ...
        'The big sparse builder requires actual per-pair Cmap data.');
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
stats.reason = 'not_requested';
stats.n_pairs = n_pairs;
stats.N_c = N_c;
stats.N_check = N_check;
stats.used_pair_cache = false;
stats.direct_source_corr = plan.direct_source_corr;
stats.direct_u_corr = plan.direct_u_corr;
stats.structured_apply = plan.structured_apply;
stats.combined_pair_u_corr = plan.combined_pair_u_corr;
stats.local_source_entries = counts.source/max(1,n_pairs);
stats.local_u_entries = counts.u/max(1,n_pairs);
stats.nnz_source = 0;
stats.nnz_u = 0;
stats.nnz_pair_nonp = 0;
stats.nnz_pair_nonp_u_corr = 0;
stats.nnz_pair_projector = 0;
stats.nnz_source_scatter = 0;
stats.nnz_u_cross = 0;
stats.nnz_u_peanut = 0;
stats.rotations_used = false;
stats.estimated_sparse_bytes = ram_estimate.estimated_sparse_bytes;
stats.estimated_auxiliary_bytes = ram_estimate.estimated_auxiliary_bytes;
stats.estimated_build_bytes = ram_estimate.estimated_build_bytes;
stats.estimated_peak_bytes = ram_estimate.estimated_peak_bytes;
stats.estimated_sparse_MB = ram_estimate.estimated_sparse_MB;
stats.estimated_auxiliary_MB = ram_estimate.estimated_auxiliary_MB;
stats.estimated_build_MB = ram_estimate.estimated_build_MB;
stats.estimated_peak_MB = ram_estimate.estimated_peak_MB;
stats.build_time = 0;
end

function big_sparse = emptyBigSparse(n_source_rows,n_source_cols,n_u_rows, ...
    n_u_cols,pair_total_rows,plan)
big_sparse = struct();
big_sparse.matrix_plan = plan;
if plan.need_dense_pair_projector
    big_sparse.P_pair = sparse(0,0);
end
if plan.need_source_scatter_indices
    big_sparse.source_scatter_rows = zeros(0,1);
end
if plan.direct_source_corr
    big_sparse.M_source_corr = sparse(n_source_rows,n_source_cols);
end
if plan.direct_u_corr
    if plan.combined_pair_u_corr
        big_sparse.M_pair_nonp_u_corr = sparse(pair_total_rows+n_u_rows, ...
            n_u_cols);
        big_sparse.pair_nonp_rows = pair_total_rows;
    else
        big_sparse.M_u_corr = sparse(n_u_rows,n_u_cols);
    end
end
if plan.need_pair_nonp_matrix
    big_sparse.M_pair_nonp = sparse(pair_total_rows,n_source_cols);
end
if plan.need_pair_projector
    big_sparse.M_pair_projector = sparse(pair_total_rows,pair_total_rows);
end
if plan.need_source_scatter
    big_sparse.M_source_scatter = sparse(n_source_rows,2*pair_total_rows);
end
if plan.need_u_cross
    big_sparse.M_u_cross = sparse(n_u_rows,n_u_cols);
end
if plan.need_u_peanut
    big_sparse.M_u_peanut = sparse(n_u_rows,pair_total_rows);
end
end

function entries = preallocateSparseEntries(counts)
entries = struct();
entries.source = newSparseEntries(counts.source);
entries.u = newSparseEntries(counts.u);
entries.pair_nonp = newSparseEntries(counts.pair_nonp);
entries.pair_nonp_u_corr = newSparseEntries(counts.pair_nonp_u_corr);
entries.pair_projector = newSparseEntries(counts.pair_projector);
entries.source_scatter = newSparseEntries(counts.source_scatter);
entries.u_cross = newSparseEntries(counts.u_cross);
entries.u_peanut = newSparseEntries(counts.u_peanut);
end

function block = newSparseEntries(n)
block = struct('rows',zeros(n,1),'cols',zeros(n,1), ...
    'vals',zeros(n,1),'next',1);
end

function entries = appendDenseBlock(entries,name,row_idx,col_idx,A)
block = entries.(name);
row_idx = row_idx(:);
col_idx = col_idx(:);
nr = numel(row_idx);
nc = numel(col_idx);
n = nr*nc;
loc = block.next:block.next+n-1;

block.rows(loc) = repmat(row_idx,nc,1);
block.cols(loc) = repelem(col_idx,nr);
block.vals(loc) = real(A(:));
block.next = block.next + n;
entries.(name) = block;
end

function entries = appendSourceScatter(entries,source_out_idx,pair_idx,pair_total_rows)
block = entries.source_scatter;
n = numel(source_out_idx);
loc = block.next:block.next+n-1;

block.rows(loc) = source_out_idx(:);
block.cols(loc) = [pair_idx, pair_total_rows + pair_idx]';
block.vals(loc) = 1;
block.next = block.next + n;
entries.source_scatter = block;
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

function S = sparseFromEntries(block,n_rows,n_cols)
if block.next ~= numel(block.vals) + 1
    error('buildMobPeanutBigSparseStokes:EntryCountMismatch', ...
        'Internal sparse entry count mismatch.');
end
S = sparse(block.rows,block.cols,block.vals,n_rows,n_cols);
end

function P_pair = buildPairProjectionMatrix(rbase)
K = getKmat2D(rbase,0);
L = K*((K'*K)\K');
P_pair = getILpair(L);
end

function [Ucross,Ecolloc] = buildActualCoarsePairDense( ...
    q,rbase_in_c,rout_base_c,pairs,row)
i = pairs(row,1);
j = pairs(row,2);

rin_pair_c = [q(i)+rbase_in_c; q(j)+rbase_in_c];
rout_pair = [q(i)+rout_base_c; q(j)+rout_base_c];

Ecolloc = stokSLPmat(rin_pair_c,rout_pair,1);
Ucross = buildCrossPairVelocityMap(Ecolloc,numel(rbase_in_c), ...
    numel(rout_base_c));
end

function Ucross = buildCrossPairVelocityMap(Epair,N_src,N_tgt)
Ucross = zeros(size(Epair));

tgt_i_x = 1:N_tgt;
tgt_j_x = N_tgt+1:2*N_tgt;
tgt_i_y = 2*N_tgt+1:3*N_tgt;
tgt_j_y = 3*N_tgt+1:4*N_tgt;

src_i_x = 1:N_src;
src_j_x = N_src+1:2*N_src;
src_i_y = 2*N_src+1:3*N_src;
src_j_y = 3*N_src+1:4*N_src;

Ucross(tgt_i_x,[src_j_x src_j_y]) = Epair(tgt_i_x,[src_j_x src_j_y]);
Ucross(tgt_j_x,[src_i_x src_i_y]) = Epair(tgt_j_x,[src_i_x src_i_y]);
Ucross(tgt_i_y,[src_j_x src_j_y]) = Epair(tgt_i_y,[src_j_x src_j_y]);
Ucross(tgt_j_y,[src_i_x src_i_y]) = Epair(tgt_j_y,[src_i_x src_i_y]);

Ucross = -Ucross;
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

function idx = pairVelocityOutputIndices(i,j,N_check,P)
pm_check = P*N_check;
idx_i = (i-1)*N_check+1:i*N_check;
idx_j = (j-1)*N_check+1:j*N_check;
idx = [idx_i, idx_j, pm_check+idx_i, pm_check+idx_j]';
end
