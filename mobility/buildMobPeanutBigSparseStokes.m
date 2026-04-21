function [big_sparse,stats] = buildMobPeanutBigSparseStokes(geom,basis)
%BUILDMOBPEANUTBIGSPARSESTOKES Build solve-grid sparse pair-correction maps.
%
% Calling buildMobPeanutBigSparseStokes with no inputs runs a visual
% self-test on a P=20 random_discs_mc packing at phi=0.65 and shows the
% sparse matrix structures with spy.
%
% This v1 builder is deliberately actual-pair only. Rotations are not part
% of this path: solve_mob_peanut_enhanced disables
% opt.reuse_pair_basis_by_sep when opt.use_big_sparse is true, so every
% close pair contributes the dense maps for its actual geometry.
%
% Fixed matrix set:
%   M_pair_nonp maps projected one-body coarse sources to unprojected
%       close-pair coarse-source increments.
%   P_pair is the dense local rigid projection applied to each pair block.
%   source_scatter_rows accumulates [projected; unprojected] pair blocks
%       into global coarse-source correction vectors.
%   M_u_corr maps projected one-body coarse sources directly to u_corr.
%       This is the default opt.big_sparse_direct_u_corr=true path.
%   M_u_cross and M_u_peanut are built instead when
%       opt.big_sparse_direct_u_corr=false for comparison.
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
rbase_in_c = geom.rbase_in_c(:);
rcheck = geom.rcheck;
N_c = opt.N_c;
P = numel(q);
n_pairs = size(pairs,1);
n_coarse = P*N_c;
N_check = numel(rcheck)/P;

validateBigSparseInputs(geom,basis,N_check);

n_source_cols = 2*n_coarse;
n_u_rows = 2*P*N_check;
n_u_cols = 2*n_coarse;
pair_block_rows = 4*N_c;
pair_total_rows = n_pairs*pair_block_rows;
P_pair = buildPairProjectionMatrix(rbase_in_c);
ram_estimate = estimateMobPeanutBigSparseRamStokes(P,N_c,N_check,n_pairs,opt);
plan = ram_estimate.matrix_plan;
counts = ram_estimate.counts;
direct_u_corr = plan.direct_u_corr;

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
    big_sparse = emptyBigSparse(n_source_cols,n_u_rows,n_u_cols, ...
        pair_total_rows,plan);
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
    u_out_idx = pairVelocityOutputIndices(i,j,N_check,P);
    pair_idx = (row-1)*pair_block_rows+1:row*pair_block_rows;

    C_nonp = basis.Cmap{i,j};
    [Ucross,Ecolloc] = buildActualCoarsePairDense( ...
        q,rbase_in_c,rout_base_c,pairs,row);

    entries = appendDenseBlock(entries,'pair_nonp',pair_idx,in_idx,C_nonp);
    if direct_u_corr
        C_proj = P_pair*C_nonp;
        entries = appendDenseBlock(entries,'u',u_out_idx,in_idx, ...
            Ucross - Ecolloc*C_proj);
    else
        entries = appendDenseBlock(entries,'u_cross',u_out_idx,in_idx, ...
            Ucross);
        entries = appendDenseBlock(entries,'u_peanut',u_out_idx,pair_idx, ...
            Ecolloc);
    end
end

big_sparse = struct();
big_sparse.matrix_plan = plan;
big_sparse.ram_estimate = ram_estimate;
big_sparse.P_pair = P_pair;
big_sparse.source_scatter_rows = buildSourceScatterRows(pairs,N_c,P);
big_sparse.M_pair_nonp = sparseFromEntries(entries.pair_nonp, ...
    pair_total_rows,n_source_cols);
if direct_u_corr
    big_sparse.M_u_corr = sparseFromEntries(entries.u,n_u_rows,n_u_cols);
else
    big_sparse.M_u_cross = sparseFromEntries(entries.u_cross,n_u_rows, ...
        n_u_cols);
    big_sparse.M_u_peanut = sparseFromEntries(entries.u_peanut,n_u_rows, ...
        pair_total_rows);
end

if direct_u_corr
    stats.nnz_u = nnz(big_sparse.M_u_corr);
else
    stats.nnz_u_cross = nnz(big_sparse.M_u_cross);
    stats.nnz_u_peanut = nnz(big_sparse.M_u_peanut);
end
stats.nnz_pair_nonp = nnz(big_sparse.M_pair_nonp);
stats.nnz_source_scatter = numel(big_sparse.source_scatter_rows);
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
stats.source_correction = 'factored_structured';
stats.velocity_correction = plan.velocity_correction;
stats.direct_u_corr = plan.direct_u_corr;
stats.local_pair_nonp_entries = counts.pair_nonp/max(1,n_pairs);
stats.local_u_entries = counts.u/max(1,n_pairs);
stats.local_u_cross_entries = counts.u_cross/max(1,n_pairs);
stats.local_u_peanut_entries = counts.u_peanut/max(1,n_pairs);
stats.nnz_u = 0;
stats.nnz_u_cross = 0;
stats.nnz_u_peanut = 0;
stats.nnz_pair_nonp = 0;
stats.nnz_source_scatter = 0;
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

function big_sparse = emptyBigSparse(n_source_cols,n_u_rows,n_u_cols, ...
    pair_total_rows,plan)
big_sparse = struct();
big_sparse.matrix_plan = plan;
big_sparse.P_pair = sparse(0,0);
big_sparse.source_scatter_rows = zeros(0,1);
big_sparse.M_pair_nonp = sparse(pair_total_rows,n_source_cols);
if plan.direct_u_corr
    big_sparse.M_u_corr = sparse(n_u_rows,n_u_cols);
else
    big_sparse.M_u_cross = sparse(n_u_rows,n_u_cols);
    big_sparse.M_u_peanut = sparse(n_u_rows,pair_total_rows);
end
end

function entries = preallocateSparseEntries(counts)
entries = struct();
entries.u = newSparseEntries(counts.u);
entries.pair_nonp = newSparseEntries(counts.pair_nonp);
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

function test_buildMobPeanutBigSparseStokes
fprintf('buildMobPeanutBigSparseStokes self-test: random_discs_mc spy plots\n');

P = 20;
phi = 0.65;
rad = 1;
geom_opt = struct('domain','boxed','phi',phi,'rad',rad, ...
    'min_gap',1e-3,'n_sweeps',200,'rng_seed',210421, ...
    'visualise',false);
[q,meta] = random_discs_mc(P,geom_opt);

opt = get2Dparams(P,24,60);
opt.delta_pair = 0.2;
opt.N_peanut = 120;
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

geom_direct = geom;
geom_direct.opt.big_sparse_direct_u_corr = true;
[big_direct,stats_direct] = buildMobPeanutBigSparseStokes( ...
    geom_direct,basis);

geom_factored = geom;
geom_factored.opt.big_sparse_direct_u_corr = false;
[big_factored,stats_factored] = buildMobPeanutBigSparseStokes( ...
    geom_factored,basis);

fprintf('  P=%d, target phi=%.3f, actual phi=%.3f, close pairs=%d\n', ...
    P,phi,meta.phi,stats_direct.n_pairs);
fprintf('  direct u:   M_pair_nonp nnz=%d, M_u_corr nnz=%d\n', ...
    stats_direct.nnz_pair_nonp,stats_direct.nnz_u);
fprintf('  factored u: M_pair_nonp nnz=%d, M_u_cross nnz=%d, M_u_peanut nnz=%d\n', ...
    stats_factored.nnz_pair_nonp,stats_factored.nnz_u_cross, ...
    stats_factored.nnz_u_peanut);

figure('Name','buildMobPeanutBigSparseStokes spy self-test','Color','w');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
spy(big_direct.M_pair_nonp);
title('direct: M_pair_nonp','Interpreter','none');

nexttile;
spy(big_direct.M_u_corr);
title('direct: M_u_corr','Interpreter','none');

nexttile;
axis off
title(sprintf('P=%d, phi=%.2f, pairs=%d',P,meta.phi,stats_direct.n_pairs));

nexttile;
spy(big_factored.M_pair_nonp);
title('factored: M_pair_nonp','Interpreter','none');

nexttile;
spy(big_factored.M_u_cross);
title('factored: M_u_cross','Interpreter','none');

nexttile;
spy(big_factored.M_u_peanut);
title('factored: M_u_peanut','Interpreter','none');
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
