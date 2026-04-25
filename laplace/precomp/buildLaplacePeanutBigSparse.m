function [big_sparse,stats,pair_cache] = buildLaplacePeanutBigSparse(...
    geom,basis)
%BUILDLAPLACEPEANUTBIGSPARSE Build solve-grid sparse peanut correction maps.
%
% Calling buildLaplacePeanutBigSparse with no inputs runs a visual self-test
% on a three-body shared-separation Laplace case and shows sparse matrix
% structures with spy.

if nargin == 0
    test_buildLaplacePeanutBigSparse;
    big_sparse = [];
    stats = [];
    pair_cache = [];
    return
end

if nargin < 2
    basis = struct();
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
project_charge = logical(getOptField(opt,'project_charge',false));
build_mode = lower(char(getOptField(opt,'lap_big_sparse_build_mode','auto')));
switch build_mode
    case 'auto'
        if logical(getOptField(opt,'get_bndry_field',true))
            build_mode = 'precomputed';
        else
            build_mode = 'streaming';
        end
    case {'precomputed','streaming'}
        % supported
    otherwise
        error('buildLaplacePeanutBigSparse:BadMode',...
            ['opt.lap_big_sparse_build_mode must be ''auto'', ',...
             '''precomputed'', or ''streaming''.']);
end

validateBigSparseInputs(geom,N_check,build_mode);

ram_estimate = estimateLaplacePeanutBigSparseRam(P,N_c,N_check,n_pairs,opt);
stats = initLaplaceBigSparseStats(n_pairs,N_c,N_check,ram_estimate,...
    build_mode);
stats.requested = true;

build_limit = resolveLaplaceBigSparseBuildMemoryLimit(opt);
max_build_bytes = build_limit.max_build_bytes;
stats.max_build_bytes = build_limit.max_build_bytes;
stats.max_build_MB = build_limit.max_build_bytes/1024^2;
stats.max_build_ram_fraction = build_limit.ram_fraction;
stats.available_ram_bytes = build_limit.available_ram_bytes;
stats.available_ram_MB = build_limit.available_ram_bytes/1024^2;
stats.guardrail_scope = build_limit.scope;
stats.guardrail_checked_bytes = stats.estimated_peak_bytes;
stats.guardrail_checked_MB = stats.estimated_peak_MB;
if stats.guardrail_checked_bytes > max_build_bytes
    error('buildLaplacePeanutBigSparse:BuildMemoryLimit',...
        ['Estimated sparse-build peak memory %.3g bytes exceeds ',...
         'Laplace big-sparse build guardrail %.3g bytes.'],...
        stats.guardrail_checked_bytes,max_build_bytes);
end

pair_rows = 2*N_c;
pair_total_rows = n_pairs*pair_rows;
P_pair = buildChargeProjectionPairMatrix(N_c,project_charge);

if n_pairs == 0
    big_sparse = emptyBigSparse(pair_total_rows,n_coarse,P*N_check,...
        2*n_pairs,ram_estimate,P_pair);
    pair_cache = emptyPairCacheForMode(build_mode,opt);
    stats.active = true;
    stats.reason = 'no_close_pairs';
    stats.build_time = toc(timer);
    big_sparse.stats = stats;
    return
end

entries = initSparseBuilders(P,N_c,N_check,n_pairs,opt);

switch build_mode
    case 'precomputed'
        [entries,pair_cache,stats] = appendPrecomputedBlocks(entries,...
            geom,basis,P_pair,stats);
    case 'streaming'
        [entries,pair_cache,stats] = appendStreamingBlocks(entries,...
            geom,P_pair,stats);
    otherwise
        error('buildLaplacePeanutBigSparse:BadMode',...
            'Unsupported build mode "%s".',build_mode);
end

entries = flushAllSparseBuilders(entries);

%% Assemble output struct
big_sparse              = struct();
big_sparse.matrix_plan  = ram_estimate.matrix_plan;
big_sparse.ram_estimate = ram_estimate;
big_sparse.P_pair       = P_pair;
big_sparse.N_c          = N_c;
big_sparse.pairs        = pairs;
big_sparse.pair_rotation_meta = pair_cache.meta;

% Phase matrix cache: eliminates per-matvec struct-field extraction in
% rotatePairInputs / rotatePairOutputs when pair_map_mode='rotated_factored'.
% Built once here; costs 2*N_c*n_pairs*16 bytes (negligible vs sparse mats).
big_sparse.phase_c_mat     = [];
big_sparse.phase_c_inv_mat = [];
if n_pairs > 0 && strcmp(entries.pair_map_mode,'rotated_factored')
    meta_all = pair_cache.meta;
    if isnumeric(meta_all(1).phase_c)
        phase_c_mat     = zeros(N_c, n_pairs);
        phase_c_inv_mat = zeros(N_c, n_pairs);
        for row = 1:n_pairs
            phase_c_mat(:,row)     = meta_all(row).phase_c(:);
            phase_c_inv_mat(:,row) = meta_all(row).phase_c_inv(:);
        end
        big_sparse.phase_c_mat     = phase_c_mat;
        big_sparse.phase_c_inv_mat = phase_c_inv_mat;
    end
end

big_sparse.source_scatter_rows = buildSourceScatterRows(pairs,N_c,P);
big_sparse.qv_scatter_rows     = buildQvScatterRows(pairs);
big_sparse.pair_map_mode       = entries.pair_map_mode;

if strcmp(entries.pair_map_mode,'rotated_factored')
    big_sparse.M_pair_nonp_canon = entries.pair_nonp_canon.S;
    big_sparse.M_pair_qv_canon   = entries.qv_canon.S;
else
    big_sparse.M_pair_nonp = entries.pair_nonp.S;
    big_sparse.M_pair_qv   = entries.qv.S;
end
big_sparse.M_u_cross  = entries.u_cross.S;
big_sparse.M_u_peanut = entries.u_peanut.S;
big_sparse.M_u_qv     = entries.u_qv.S;

if strcmp(big_sparse.pair_map_mode,'rotated_factored')
    stats.nnz_pair_nonp       = 0;
    stats.nnz_pair_nonp_canon = nnz(big_sparse.M_pair_nonp_canon);
    stats.nnz_qv              = 0;
    stats.nnz_qv_canon        = nnz(big_sparse.M_pair_qv_canon);
else
    stats.nnz_pair_nonp       = nnz(big_sparse.M_pair_nonp);
    stats.nnz_pair_nonp_canon = 0;
    stats.nnz_qv              = nnz(big_sparse.M_pair_qv);
    stats.nnz_qv_canon        = 0;
end
stats.nnz_u_cross        = nnz(big_sparse.M_u_cross);
stats.nnz_u_peanut       = nnz(big_sparse.M_u_peanut);
stats.nnz_u_qv           = nnz(big_sparse.M_u_qv);
stats.nnz_u              = stats.nnz_u_cross + stats.nnz_u_peanut + stats.nnz_u_qv;
stats.nnz_source_scatter = numel(big_sparse.source_scatter_rows);
stats.nnz_qv_scatter     = numel(big_sparse.qv_scatter_rows);
stats.active             = true;
stats.reason             = '';
stats.build_time         = toc(timer);
big_sparse.stats         = stats;
end

% =========================================================================
%  LOCAL FUNCTIONS  (unchanged from original except where noted)
% =========================================================================

function validateBigSparseInputs(geom,N_check,build_mode)
opt = geom.opt;
if ~logical(getOptField(opt,'cmap',false))
    error('buildLaplacePeanutBigSparse:UnsupportedOption',...
        'opt.use_big_sparse=1 requires opt.cmap=1.');
end
if getOptField(opt,'N_peanut',0) <= 0
    error('buildLaplacePeanutBigSparse:UnsupportedOption',...
        'opt.use_big_sparse=1 requires opt.N_peanut > 0.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('buildLaplacePeanutBigSparse:UnsupportedGrid',...
        'The Laplace big-sparse matvec is only supported on the solve grid.');
end
if N_check ~= round(N_check)
    error('buildLaplacePeanutBigSparse:BadGridSize',...
        'rcheck length must be divisible by the number of particles.');
end
if strcmp(build_mode,'streaming') &&...
        logical(getOptField(opt,'get_bndry_field',false))
    error('buildLaplacePeanutBigSparse:StreamingBoundaryUnsupported',...
        ['lap_big_sparse_build_mode=''streaming'' is only supported when ',...
         'opt.get_bndry_field=0. Use ''auto'' or ''precomputed'' for ',...
         'boundary postprocessing.']);
end
if ~isfield(geom,'rout_base_f') || isempty(geom.rout_base_f)
    error('buildLaplacePeanutBigSparse:MissingFineOutputGrid',...
        'geom.rout_base_f is required for Laplace big-sparse builds.');
end
end

function stats = initLaplaceBigSparseStats(n_pairs,N_c,N_check,...
    ram_estimate,build_mode)
counts = ram_estimate.counts;
stats = struct();
stats.requested = false;
stats.active = false;
stats.backend = 'global_block_sparse';
stats.reason = 'not_requested';
stats.build_mode = build_mode;
stats.n_pairs = n_pairs;
stats.N_c = N_c;
stats.N_check = N_check;
stats.used_pair_cache = false;
stats.rotations_used = false;
stats.source_correction = 'factored_structured';
stats.field_correction = 'factored_sparse';
stats.local_pair_nonp_entries = counts.pair_nonp/max(1,n_pairs);
stats.local_u_entries = counts.u/max(1,n_pairs);
stats.local_qv_entries = counts.qv/max(1,n_pairs);
stats.nnz_u = 0;
stats.nnz_u_cross = 0;
stats.nnz_u_peanut = 0;
stats.nnz_u_qv = 0;
stats.nnz_qv = 0;
stats.nnz_qv_canon = 0;
stats.nnz_pair_nonp = 0;
stats.nnz_pair_nonp_canon = 0;
stats.nnz_source_scatter = 0;
stats.nnz_qv_scatter = 0;
stats.estimated_sparse_bytes = ram_estimate.estimated_sparse_bytes;
stats.estimated_auxiliary_bytes = ram_estimate.estimated_auxiliary_bytes;
stats.estimated_build_bytes = ram_estimate.estimated_build_bytes;
stats.estimated_peak_bytes = ram_estimate.estimated_peak_bytes;
stats.estimated_sparse_MB = ram_estimate.estimated_sparse_MB;
stats.estimated_auxiliary_MB = ram_estimate.estimated_auxiliary_MB;
stats.estimated_build_MB = ram_estimate.estimated_build_MB;
stats.estimated_peak_MB = ram_estimate.estimated_peak_MB;
stats.max_build_bytes = inf;
stats.max_build_MB = inf;
stats.max_build_ram_fraction = [];
stats.available_ram_bytes = nan;
stats.available_ram_MB = nan;
stats.guardrail_scope = 'estimated_peak_bytes';
stats.guardrail_checked_bytes = ram_estimate.estimated_peak_bytes;
stats.guardrail_checked_MB = ram_estimate.estimated_peak_MB;
stats.build_time = 0;
end

function build_limit = resolveLaplaceBigSparseBuildMemoryLimit(opt)
explicit_bytes = getOptField(opt,'lap_big_sparse_max_build_bytes',...
    getOptField(opt,'big_sparse_max_build_bytes',inf));
ram_fraction = getOptField(opt,'lap_big_sparse_max_build_ram_fraction',[]);
if isempty(explicit_bytes)
    explicit_bytes = inf;
end
if isempty(ram_fraction) || (isnumeric(ram_fraction) && isnan(ram_fraction))
    available_bytes = nan;
    max_build_bytes = explicit_bytes;
else
    validateattributes(ram_fraction,{'numeric'},...
        {'scalar','real','positive','<=',1},mfilename,...
        'opt.lap_big_sparse_max_build_ram_fraction');
    available_bytes = getLinuxMemAvailableBytes();
    max_build_bytes = min(explicit_bytes,ram_fraction*available_bytes);
end
build_limit = struct('max_build_bytes',max_build_bytes,...
    'ram_fraction',ram_fraction,...
    'available_ram_bytes',available_bytes,...
    'scope','estimated_peak_bytes');
end

function bytes = getLinuxMemAvailableBytes()
meminfo_file = '/proc/meminfo';
if exist(meminfo_file,'file') ~= 2
    error('buildLaplacePeanutBigSparse:MemAvailableUnsupported',...
        ['opt.lap_big_sparse_max_build_ram_fraction requires ',...
         '/proc/meminfo. Set opt.lap_big_sparse_max_build_bytes directly ',...
         'on this platform. The guardrail is checked against ',...
         'estimated_peak_bytes.']);
end
txt = fileread(meminfo_file);
tok = regexp(txt,'MemAvailable:\s+(\d+)\s+kB','tokens','once');
if isempty(tok)
    error('buildLaplacePeanutBigSparse:MemAvailableMissing',...
        'Could not read MemAvailable from /proc/meminfo.');
end
bytes = str2double(tok{1})*1024;
end

function big_sparse = emptyBigSparse(pair_total_rows,n_coarse,n_u_rows,...
    n_qv_rows,ram_estimate,P_pair)
big_sparse = struct();
big_sparse.matrix_plan = ram_estimate.matrix_plan;
big_sparse.ram_estimate = ram_estimate;
big_sparse.P_pair = P_pair;
big_sparse.phase_c_mat = [];
big_sparse.phase_c_inv_mat = [];
big_sparse.source_scatter_rows = zeros(0,1);
big_sparse.qv_scatter_rows = zeros(0,1);
big_sparse.pair_map_mode = 'direct';
big_sparse.M_pair_nonp = sparse(pair_total_rows,n_coarse);
big_sparse.M_u_cross = sparse(n_u_rows,n_coarse);
big_sparse.M_u_peanut = sparse(n_u_rows,pair_total_rows);
big_sparse.M_u_qv = sparse(n_u_rows,n_qv_rows);
big_sparse.M_pair_qv = sparse(n_qv_rows,n_coarse);
end

function pair_cache = emptyPairCacheForMode(build_mode,opt)
pair_cache = initLaplacePairCache();
pair_cache.enabled = logical(getOptField(opt,'reuse_pair_basis_by_sep',false));
pair_cache.stats.branch = ternary(pair_cache.enabled,'canonical_group','per_pair');
pair_cache.stats.n_pairs = 0;
pair_cache.stats.n_groups = 0;
pair_cache.stats.used_parallel = false;
pair_cache.stats.build_mode = build_mode;
end

function entries = initSparseBuilders(P,N_c,N_check,n_pairs,opt)
chunk_pairs = max(1,round(getOptField(opt,'lap_big_sparse_chunk_pairs',8)));
chunk_pairs = min(max(1,n_pairs),chunk_pairs);
n_coarse = P*N_c;
pair_rows = 2*N_c;
use_factored_pair_maps = logical(getOptField(opt,...
    'reuse_pair_basis_by_sep',false));

entries = struct();
entries.pair_map_mode = ternary(use_factored_pair_maps,'rotated_factored',...
    'direct');
if use_factored_pair_maps
    entries.pair_nonp = newSparseBuilder(0,n_pairs*pair_rows,n_coarse);
    entries.pair_nonp_canon = newSparseBuilder(...
        chunk_pairs*pair_rows*pair_rows,n_pairs*pair_rows,...
        n_pairs*pair_rows);
else
    entries.pair_nonp = newSparseBuilder(...
        chunk_pairs*pair_rows*(2*N_c),n_pairs*pair_rows,n_coarse);
    entries.pair_nonp_canon = newSparseBuilder(0,n_pairs*pair_rows,...
        n_pairs*pair_rows);
end
entries.u_cross = newSparseBuilder(...
    chunk_pairs*(2*N_check)*(2*N_c),P*N_check,n_coarse);
entries.u_peanut = newSparseBuilder(...
    chunk_pairs*(2*N_check)*pair_rows,P*N_check,n_pairs*pair_rows);
entries.u_qv = newSparseBuilder(...
    chunk_pairs*(2*N_check)*2,P*N_check,2*n_pairs);
if use_factored_pair_maps
    entries.qv = newSparseBuilder(0,2*n_pairs,n_coarse);
    entries.qv_canon = newSparseBuilder(...
        chunk_pairs*2*pair_rows,2*n_pairs,n_pairs*pair_rows);
else
    entries.qv = newSparseBuilder(...
        chunk_pairs*2*(2*N_c),2*n_pairs,n_coarse);
    entries.qv_canon = newSparseBuilder(0,2*n_pairs,n_pairs*pair_rows);
end
end

function block = newSparseBuilder(capacity,n_rows,n_cols)
capacity = max(0,round(capacity));
block = struct();
block.rows = zeros(capacity,1);
block.cols = zeros(capacity,1);
block.vals = zeros(capacity,1);
block.next = 1;
block.S = sparse(n_rows,n_cols);
end

function [entries,pair_cache,stats] = appendPrecomputedBlocks(entries,...
    geom,basis,P_pair,stats)
pairs = geom.pairs;
q = geom.q(:);
opt = geom.opt;
N_c = opt.N_c;
n_pairs = size(pairs,1);
N_check = numel(geom.rcheck)/numel(q);
rout_base_c = geom.rcheck(1:N_check)-q(1);
project_charge = logical(getOptField(opt,'project_charge',false));

if ~isfield(basis,'pair_cache') || isempty(basis.pair_cache)
    error('buildLaplacePeanutBigSparse:MissingPairCache',...
        'precomputed mode requires basis.pair_cache.');
end
pair_cache = basis.pair_cache;
use_pair_cache = isfield(pair_cache,'enabled') && pair_cache.enabled;
stats.used_pair_cache = use_pair_cache;
stats.rotations_used = use_pair_cache;

for row = 1:n_pairs
    i = pairs(row,1);
    j = pairs(row,2);
    meta = pair_cache.meta(row);
    if use_pair_cache
        group = pair_cache.groups(meta.group_id);
        [Ucross,Ec,Lr] = getActualCollocFactors(meta,q,geom.rbase_in_c,...
            rout_base_c);
        entries = appendFactoredCanonicalPairBlocks(entries,pairs,row,...
            N_c,N_check,numel(q),group,Ucross,Ec,Lr,P_pair,...
            project_charge);
    else
        C_nonp = basis.Cmap{i,j};
        QV = basis.Cmap_QV{i,j};
        [Ucross,Ec,Lr] = getActualCollocFactors(meta,q,geom.rbase_in_c,...
            rout_base_c);
        entries = appendPairBlocks(entries,pairs,row,N_c,N_check,numel(q),...
            C_nonp,QV,Ucross,Ec,Lr,P_pair,project_charge);
    end
end
pair_cache.stats.build_mode = stats.build_mode;
end

function [entries,pair_cache,stats] = appendStreamingBlocks(entries,...
    geom,P_pair,stats)
q = geom.q(:);
pairs = geom.pairs;
opt = geom.opt;
N_c = opt.N_c;
N_f = opt.N_f;
n_pairs = size(pairs,1);
N_check = numel(geom.rcheck)/numel(q);
rout_base_c = geom.rcheck(1:N_check)-q(1);
project_charge = logical(getOptField(opt,'project_charge',false));
reuse_pair_basis = logical(getOptField(opt,'reuse_pair_basis_by_sep',false));
show_counter = logical(getOptField(opt,'show_counter',false));

pair_cache = initLaplacePairCache();
pair_cache.enabled = reuse_pair_basis;
pair_cache.check_rotations = false;
pair_cache.shared_sep_tol = getOptField(opt,'shared_sep_tol',0);
pair_cache.rout_base_f = geom.rout_base_f(:);
pair_cache.meta = buildLaplacePairMeta(q,pairs,N_c,N_f,opt);
pair_cache.stats.n_pairs = n_pairs;
pair_cache.stats.used_parallel = false;
pair_cache.stats.build_mode = stats.build_mode;
stats.used_pair_cache = reuse_pair_basis;
stats.rotations_used = reuse_pair_basis;

if ~reuse_pair_basis
    pair_cache.stats.branch = 'per_pair';
    for row = 1:n_pairs
        group = buildLaplacePairGroup([],row,pair_cache.meta(row).sep,...
            q,geom.rbase_in_c,geom.rbase_in_f,geom.rimage_vec,...
            geom.refine,pairs,opt,geom.rout_base_f,project_charge,false);
        [Ucross,Ec,Lr] = buildLaplaceActualPairCollocFactors(...
            pair_cache.meta(row),q,geom.rbase_in_c,rout_base_c);
        entries = appendPairBlocks(entries,pairs,row,N_c,N_check,numel(q),...
            group.Cmap,group.Cmap_QV,Ucross,Ec,Lr,P_pair,...
            project_charge);
        if show_counter
            fprintf('buildLaplacePeanutBigSparse: streamed pair %d/%d\n',...
                row,n_pairs);
        end
    end
    return
end

[group_id,group_sep,rep_rows] = groupLaplacePairSeparations(...
    pair_cache.meta,pair_cache.shared_sep_tol);
pair_cache.group_id = group_id;
pair_cache.group_sep = group_sep;
pair_cache.representative_rows = rep_rows;
pair_cache.n_groups = numel(group_sep);
pair_cache.stats.branch = 'canonical_group';
pair_cache.stats.n_groups = pair_cache.n_groups;

for row = 1:n_pairs
    gid = group_id(row);
    pair_cache.meta(row).group_id = gid;
    pair_cache.meta(row).sep = group_sep(gid);
end

for gg = 1:pair_cache.n_groups
    group = buildLaplacePairGroup(gg,rep_rows(gg),group_sep(gg),...
        q,geom.rbase_in_c,geom.rbase_in_f,geom.rimage_vec,...
        geom.refine,pairs,opt,geom.rout_base_f,project_charge,true);
    rows = find(group_id == gg).';
    block_cache = containers.Map('KeyType','char','ValueType','any');
    for row = rows
        meta = pair_cache.meta(row);
        cache_key = pairBlockCacheKey(meta.rot);
        if isKey(block_cache,cache_key)
            block = block_cache(cache_key);
            Ucross = block.Ucross;
            Ec = block.Ec;
            Lr = block.Lr;
        else
            [Ucross,Ec,Lr] = buildLaplaceActualPairCollocFactors(...
                meta,q,geom.rbase_in_c,rout_base_c);
            block_cache(cache_key) = struct('Ucross',Ucross,'Ec',Ec,...
                'Lr',Lr);
        end
        entries = appendFactoredCanonicalPairBlocks(entries,pairs,row,N_c,...
            N_check,numel(q),group,Ucross,Ec,Lr,P_pair,project_charge);
    end
    if show_counter
        fprintf(['buildLaplacePeanutBigSparse: streamed canonical group ',...
            '%d/%d covering %d pairs\n'],gg,pair_cache.n_groups,numel(rows));
    end
end
end

function key = pairBlockCacheKey(rot)
key = sprintf('%.16g_%.16g',real(rot),imag(rot));
end

function [Ucross,Ec,Lr] = getActualCollocFactors(meta,q,rbase_in_c,...
    rout_base_c)
if isfield(meta,'Ucross_colloc_actual') &&...
        ~isempty(meta.Ucross_colloc_actual) &&...
        isfield(meta,'Ec_colloc_actual') && ~isempty(meta.Ec_colloc_actual)
    Ucross = meta.Ucross_colloc_actual;
    Ec = meta.Ec_colloc_actual;
    if isfield(meta,'Lr_colloc_actual') && ~isempty(meta.Lr_colloc_actual)
        Lr = meta.Lr_colloc_actual;
    else
        Lr = buildLr(numel(rout_base_c));
    end
else
    [Ucross,Ec,Lr] = buildLaplaceActualPairCollocFactors(...
        meta,q,rbase_in_c,rout_base_c);
end
end

function Lr = buildLr(n)
Lr = [ones(n,1) zeros(n,1); zeros(n,1) ones(n,1)];
end

function entries = appendPairBlocks(entries,pairs,row,N_c,N_check,P,...
    C_nonp,QV,Ucross,Ec,Lr,P_pair,project_charge)
i = pairs(row,1);
j = pairs(row,2);
in_idx = pairCoarseInputIndices(i,j,N_c);
global_in_idx = pairCoarseInputIndicesGlobal(i,j,N_c,P);
pair_idx = (row-1)*(2*N_c)+1:row*(2*N_c);
u_idx = pairOutputIndices(i,j,N_check);
qv_idx = (row-1)*2+1:row*2;

entries.pair_nonp = appendDenseBlock(entries.pair_nonp,pair_idx,...
    global_in_idx,C_nonp);
entries.u_cross = appendDenseBlock(entries.u_cross,u_idx,global_in_idx,...
    Ucross);
entries.u_peanut = appendDenseBlock(entries.u_peanut,u_idx,pair_idx,Ec);
if project_charge
    entries.u_qv = appendDenseBlock(entries.u_qv,u_idx,qv_idx,Lr);
end
entries.qv = appendDenseBlock(entries.qv,qv_idx,global_in_idx,QV);

if numel(in_idx) ~= 2*N_c
    error('buildLaplacePeanutBigSparse:BadLocalIndex',...
        'Internal pair input index size mismatch.');
end
end

function entries = appendFactoredCanonicalPairBlocks(entries,pairs,row,N_c,...
    N_check,P,group,Ucross,Ec,Lr,P_pair,project_charge)
i = pairs(row,1);
j = pairs(row,2);
global_in_idx = pairCoarseInputIndicesGlobal(i,j,N_c,P);
pair_idx = (row-1)*(2*N_c)+1:row*(2*N_c);
u_idx = pairOutputIndices(i,j,N_check);
qv_idx = (row-1)*2+1:row*2;

entries.pair_nonp_canon = appendDenseBlock(entries.pair_nonp_canon,...
    pair_idx,pair_idx,group.Cmap);
entries.u_cross = appendDenseBlock(entries.u_cross,u_idx,global_in_idx,...
    Ucross);
entries.u_peanut = appendDenseBlock(entries.u_peanut,u_idx,pair_idx,Ec);
if project_charge
    entries.u_qv = appendDenseBlock(entries.u_qv,u_idx,qv_idx,Lr);
end
entries.qv_canon = appendDenseBlock(entries.qv_canon,qv_idx,pair_idx,...
    group.Cmap_QV);

if size(P_pair,1) ~= 2*N_c
    error('buildLaplacePeanutBigSparse:BadPairProjector',...
        'Internal pair projector size mismatch.');
end
end

function block = appendDenseBlock(block,row_idx,col_idx,A)
row_idx = row_idx(:);
col_idx = col_idx(:);
nr = numel(row_idx);
nc = numel(col_idx);
n = nr*nc;
if n == 0
    return
end
if block.next+n-1 > numel(block.vals)
    block = flushSparseBuilder(block);
end
if n > numel(block.vals)
    block.S = block.S + sparse(repmat(row_idx,nc,1),...
        repelem(col_idx,nr),real(A(:)),size(block.S,1),size(block.S,2));
    return
end
loc = block.next:block.next+n-1;
block.rows(loc) = repmat(row_idx,nc,1);
block.cols(loc) = repelem(col_idx,nr);
block.vals(loc) = real(A(:));
block.next = block.next + n;
end

function block = flushSparseBuilder(block)
n = block.next - 1;
if n > 0
    block.S = block.S + sparse(block.rows(1:n),block.cols(1:n),...
        block.vals(1:n),size(block.S,1),size(block.S,2));
    block.next = 1;
end
end

function entries = flushAllSparseBuilders(entries)
entries.pair_nonp        = flushSparseBuilder(entries.pair_nonp);
entries.pair_nonp_canon  = flushSparseBuilder(entries.pair_nonp_canon);
entries.u_cross          = flushSparseBuilder(entries.u_cross);
entries.u_peanut         = flushSparseBuilder(entries.u_peanut);
entries.u_qv             = flushSparseBuilder(entries.u_qv);
entries.qv               = flushSparseBuilder(entries.qv);
entries.qv_canon         = flushSparseBuilder(entries.qv_canon);
end

function P_pair = buildChargeProjectionPairMatrix(N_c,project_charge)
if project_charge
    P_body = eye(N_c) - ones(N_c,N_c)/N_c;
else
    P_body = eye(N_c);
end
Z = zeros(N_c,N_c);
P_pair = [P_body Z; Z P_body];
end

function rows = buildSourceScatterRows(pairs,N_c,P)
n_pairs  = size(pairs,1);
n_coarse = P*N_c;
local    = (0:N_c-1)';
i_starts = (pairs(:,1)-1)*N_c;
j_starts = (pairs(:,2)-1)*N_c;
ij_local = [bsxfun(@plus,i_starts',local);...
            bsxfun(@plus,j_starts',local)];   % (2*N_c) x n_pairs
block_offsets = [0; n_coarse];
rows = bsxfun(@plus,reshape(block_offsets,1,1,2),ij_local) + 1;
rows = reshape(permute(rows,[1 3 2]),4*N_c*n_pairs,1);
end

function rows = buildQvScatterRows(pairs)
n_pairs = size(pairs,1);
rows = zeros(2*n_pairs,1);
for row = 1:n_pairs
    rows((row-1)*2+1:row*2) = pairs(row,:).';
end
end

function idx = pairCoarseInputIndicesGlobal(i,j,N_c,P)
idx_i = (i-1)*N_c+1:i*N_c;
idx_j = (j-1)*N_c+1:j*N_c;
idx = [idx_i idx_j].';
if max(idx) > P*N_c
    error('buildLaplacePeanutBigSparse:BadGlobalIndex',...
        'Internal global input index out of range.');
end
end

function idx = pairCoarseInputIndices(i,j,N_c) %#ok<INUSD>
idx = (1:2*N_c).';
end

function idx = pairOutputIndices(i,j,N_check)
idx_i = (i-1)*N_check+1:i*N_check;
idx_j = (j-1)*N_check+1:j*N_check;
idx = [idx_i idx_j].';
end

function out = ternary(cond,a,b)
if cond; out = a; else; out = b; end
end

% =========================================================================
%  SELF-TEST  (unchanged)
% =========================================================================

function test_buildLaplacePeanutBigSparse
fprintf('buildLaplacePeanutBigSparse self-test: shared-pair spy plots\n');
rng(220422);

[geom,basis] = buildSelfTestLaplaceData();
tau = randn(numel(geom.rvec_out),1);
res_ref = matvec_lap_peanut_enhanced(tau,geom,basis);

geom_pre = geom;
geom_pre.opt.use_big_sparse = true;
geom_pre.opt.lap_big_sparse_build_mode = 'precomputed';
[big_pre,stats_pre,pair_cache_pre] = buildLaplacePeanutBigSparse(...
    geom_pre,basis);
basis_pre = basis;
basis_pre.big_sparse = big_pre;
basis_pre.pair_cache = pair_cache_pre;
geom_pre.pair_cache = pair_cache_pre;
res_pre = matvec_lap_peanut_big_sparse(tau,geom_pre,basis_pre);

geom_stream = geom_pre;
geom_stream.opt.lap_big_sparse_build_mode = 'streaming';
basis_stream = basis;
basis_stream.Upf = [];
basis_stream.Ypf = [];
basis_stream.DC_all = [];
basis_stream.YC_all = [];
basis_stream.Cmap = [];
basis_stream.Cmap_QV = [];
basis_stream.pair_cache = initLaplacePairCache();
[big_stream,stats_stream,pair_cache_stream] = buildLaplacePeanutBigSparse(...
    geom_stream,basis_stream);
basis_stream.big_sparse = big_stream;
basis_stream.pair_cache = pair_cache_stream;
geom_stream.pair_cache = pair_cache_stream;
res_stream = matvec_lap_peanut_big_sparse(tau,geom_stream,basis_stream);

err_pre    = relerrSelfTest(res_pre,res_ref);
err_stream = relerrSelfTest(res_stream,res_ref);

fprintf('  P=%d, close pairs=%d, shared groups=%d\n',...
    numel(geom.q),stats_pre.n_pairs,pair_cache_pre.n_groups);
fprintf(['  precomputed: build %.3fs, sparse %.2f MB, ',...
    'canonical pair nnz=%d, relerr %.3e\n'],...
    stats_pre.build_time,stats_pre.estimated_sparse_MB,...
    stats_pre.nnz_pair_nonp_canon,err_pre);
fprintf(['  streaming:   build %.3fs, sparse %.2f MB, ',...
    'canonical pair nnz=%d, relerr %.3e\n'],...
    stats_stream.build_time,stats_stream.estimated_sparse_MB,...
    stats_stream.nnz_pair_nonp_canon,err_stream);

figure('Name','buildLaplacePeanutBigSparse spy self-test','Color','w');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile; spy(big_pre.M_pair_nonp_canon);
title('precomputed: M_pair_nonp_canon','Interpreter','none');
nexttile; spy(big_pre.M_u_cross);
title('precomputed: M_u_cross','Interpreter','none');
nexttile; spy(big_pre.M_u_peanut);
title('precomputed: M_u_peanut','Interpreter','none');
nexttile; spy(big_stream.M_pair_nonp_canon);
title('streaming: M_pair_nonp_canon','Interpreter','none');
nexttile; spy(big_stream.M_u_cross);
title('streaming: M_u_cross','Interpreter','none');
nexttile; spy(big_stream.M_u_peanut);
title('streaming: M_u_peanut','Interpreter','none');

    function [geom,basis] = buildSelfTestLaplaceData()
    R = 2; gap = 2e-2;
    P = 50; theta = 1.3;
    q = exp(1i*theta)*grow_cluster(P,gap,2,2);
    N_c = 20; N_f = 28;
    opt = getLaplace2Dparams(P,R,N_c,N_f);
    opt.N_peanut = 64;
    opt.delta_pair = 0.2*R;
    opt.cmap = 1;
    opt.compress_cmap = 0;
    opt.reuse_pair_basis_by_sep = true;
    opt.shared_sep_tol = 1e-8;
    opt.get_bndry_field = 0;
    opt.use_fmm = false;
    opt.show_counter = 0;
    opt.visualise_sol = 0;
    opt.visualise_grid = 0;
    opt.gmres_verbose = 0;
    [rbase_in_c,rbase_out_c,rbase_in_f,rout_base_f,nout] =...
        buildSelfTestCircleGrids(R,opt);
    rvec_in_c = zeros(P*N_c,1);
    rout = zeros(P*nout,1);
    for k = 1:P
        rvec_in_c((k-1)*N_c+1:k*N_c) = q(k)+rbase_in_c;
        rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
    end
    [~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
    [UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_QV,pair_cache] =...
        getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f,...
        rbase_out_c,rimage_vec,refine,pairs,opt);
    [U,Y] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout]);
    geom = struct();
    geom.rbase_in_c  = rbase_in_c;
    geom.rbase_in_f  = rbase_in_f;
    geom.rout_base_f = rout_base_f;
    geom.refine      = refine;
    geom.opt         = opt;
    geom.rvec_out    = rout;
    geom.rcheck      = rout;
    geom.q           = q;
    geom.pairs       = pairs;
    geom.rimage_vec  = rimage_vec;
    geom.rvec_in     = rvec_in_c;
    geom.pair_cache  = pair_cache;
    basis = struct();
    basis.U          = U;
    basis.Y          = Y;
    basis.Upf        = UB_all;
    basis.Ypf        = YB_all;
    basis.DC_all     = UC_all;
    basis.YC_all     = YC_all;
    basis.Cmap       = Cmap;
    basis.Cmap_QV    = Cmap_QV;
    basis.pair_cache = pair_cache;
    basis.Nii        = lapSLPmat(rbase_in_c,rbase_out_c);
    end

    function [rbase_in_c,rbase_out_c,rbase_in_f,rout_base_f,nout] =...
            buildSelfTestCircleGrids(R,opt)
    N_c = opt.N_c; N_f = opt.N_f;
    nout = ceil(opt.a_c*N_c);
    t = linspace(0,2*pi,N_c+1)'; t = t(1:end-1);
    rbase_in_c = opt.Rp_c*exp(1i*t);
    t = linspace(0,2*pi,nout+1)'; t = t(1:end-1);
    rbase_out_c = R*exp(1i*t);
    t = linspace(0,2*pi,N_f+1)'; t = t(1:end-1);
    rbase_in_f = opt.Rp_f*exp(1i*t);
    t = linspace(0,2*pi,ceil(opt.a_f*N_f)+1)'; t = t(1:end-1);
    rout_base_f = R*exp(1i*t);
    end

    function e = relerrSelfTest(a,b)
    e = norm(a-b,inf)/max(1,norm(b,inf));
    end
end

% function [big_sparse,stats,pair_cache] = buildLaplacePeanutBigSparse( ...
%     geom,basis)
% %BUILDLAPLACEPEANUTBIGSPARSE Build solve-grid sparse peanut correction maps.
% %
% % Calling buildLaplacePeanutBigSparse with no inputs runs a visual self-test
% % on a three-body shared-separation Laplace case and shows sparse matrix
% % structures with spy.
% 
% if nargin == 0
%     test_buildLaplacePeanutBigSparse;
%     big_sparse = [];
%     stats = [];
%     pair_cache = [];
%     return
% end
% 
% if nargin < 2
%     basis = struct();
% end
% 
% timer = tic;
% 
% opt = geom.opt;
% q = geom.q(:);
% pairs = geom.pairs;
% rcheck = geom.rcheck;
% N_c = opt.N_c;
% P = numel(q);
% n_pairs = size(pairs,1);
% n_coarse = P*N_c;
% N_check = numel(rcheck)/P;
% project_charge = logical(getOptField(opt,'project_charge',false));
% build_mode = lower(char(getOptField(opt,'lap_big_sparse_build_mode','auto')));
% switch build_mode
%     case 'auto'
%         if logical(getOptField(opt,'get_bndry_field',true))
%             build_mode = 'precomputed';
%         else
%             build_mode = 'streaming';
%         end
%     case {'precomputed','streaming'}
%         % supported
%     otherwise
%         error('buildLaplacePeanutBigSparse:BadMode', ...
%             ['opt.lap_big_sparse_build_mode must be ''auto'', ', ...
%              '''precomputed'', or ''streaming''.']);
% end
% 
% validateBigSparseInputs(geom,N_check,build_mode);
% 
% ram_estimate = estimateLaplacePeanutBigSparseRam(P,N_c,N_check,n_pairs,opt);
% stats = initLaplaceBigSparseStats(n_pairs,N_c,N_check,ram_estimate, ...
%     build_mode);
% stats.requested = true;
% 
% build_limit = resolveLaplaceBigSparseBuildMemoryLimit(opt);
% max_build_bytes = build_limit.max_build_bytes;
% stats.max_build_bytes = build_limit.max_build_bytes;
% stats.max_build_MB = build_limit.max_build_bytes/1024^2;
% stats.max_build_ram_fraction = build_limit.ram_fraction;
% stats.available_ram_bytes = build_limit.available_ram_bytes;
% stats.available_ram_MB = build_limit.available_ram_bytes/1024^2;
% stats.guardrail_scope = build_limit.scope;
% stats.guardrail_checked_bytes = stats.estimated_peak_bytes;
% stats.guardrail_checked_MB = stats.estimated_peak_MB;
% if stats.guardrail_checked_bytes > max_build_bytes
%     error('buildLaplacePeanutBigSparse:BuildMemoryLimit', ...
%         ['Estimated sparse-build peak memory %.3g bytes exceeds ', ...
%          'Laplace big-sparse build guardrail %.3g bytes.'], ...
%         stats.guardrail_checked_bytes,max_build_bytes);
% end
% 
% pair_rows = 2*N_c;
% pair_total_rows = n_pairs*pair_rows;
% P_pair = buildChargeProjectionPairMatrix(N_c,project_charge);
% 
% if n_pairs == 0
%     big_sparse = emptyBigSparse(pair_total_rows,n_coarse,P*N_check, ...
%         2*n_pairs,ram_estimate,P_pair);
%     pair_cache = emptyPairCacheForMode(build_mode,opt);
%     stats.active = true;
%     stats.reason = 'no_close_pairs';
%     stats.build_time = toc(timer);
%     big_sparse.stats = stats;
%     return
% end
% 
% entries = initSparseBuilders(P,N_c,N_check,n_pairs,opt);
% 
% switch build_mode
%     case 'precomputed'
%         [entries,pair_cache,stats] = appendPrecomputedBlocks(entries, ...
%             geom,basis,P_pair,stats);
%     case 'streaming'
%         [entries,pair_cache,stats] = appendStreamingBlocks(entries, ...
%             geom,P_pair,stats);
%     otherwise
%         error('buildLaplacePeanutBigSparse:BadMode', ...
%             'Unsupported build mode "%s".',build_mode);
% end
% 
% entries = flushAllSparseBuilders(entries);
% 
% big_sparse = struct();
% big_sparse.matrix_plan = ram_estimate.matrix_plan;
% big_sparse.ram_estimate = ram_estimate;
% big_sparse.P_pair = P_pair;
% big_sparse.N_c = N_c;
% big_sparse.pairs = pairs;
% big_sparse.pair_rotation_meta = pair_cache.meta;
% big_sparse.source_scatter_rows = buildSourceScatterRows(pairs,N_c,P);
% big_sparse.qv_scatter_rows = buildQvScatterRows(pairs);
% big_sparse.pair_map_mode = entries.pair_map_mode;
% if strcmp(entries.pair_map_mode,'rotated_factored')
%     big_sparse.M_pair_nonp_canon = entries.pair_nonp_canon.S;
%     big_sparse.M_pair_qv_canon = entries.qv_canon.S;
% else
%     big_sparse.M_pair_nonp = entries.pair_nonp.S;
%     big_sparse.M_pair_qv = entries.qv.S;
% end
% big_sparse.M_u_cross = entries.u_cross.S;
% big_sparse.M_u_peanut = entries.u_peanut.S;
% big_sparse.M_u_qv = entries.u_qv.S;
% 
% if strcmp(big_sparse.pair_map_mode,'rotated_factored')
%     stats.nnz_pair_nonp = 0;
%     stats.nnz_pair_nonp_canon = nnz(big_sparse.M_pair_nonp_canon);
%     stats.nnz_qv = 0;
%     stats.nnz_qv_canon = nnz(big_sparse.M_pair_qv_canon);
% else
%     stats.nnz_pair_nonp = nnz(big_sparse.M_pair_nonp);
%     stats.nnz_pair_nonp_canon = 0;
%     stats.nnz_qv = nnz(big_sparse.M_pair_qv);
%     stats.nnz_qv_canon = 0;
% end
% stats.nnz_u_cross = nnz(big_sparse.M_u_cross);
% stats.nnz_u_peanut = nnz(big_sparse.M_u_peanut);
% stats.nnz_u_qv = nnz(big_sparse.M_u_qv);
% stats.nnz_u = stats.nnz_u_cross + stats.nnz_u_peanut + stats.nnz_u_qv;
% stats.nnz_source_scatter = numel(big_sparse.source_scatter_rows);
% stats.nnz_qv_scatter = numel(big_sparse.qv_scatter_rows);
% stats.active = true;
% stats.reason = '';
% stats.build_time = toc(timer);
% big_sparse.stats = stats;
% end
% 
% function validateBigSparseInputs(geom,N_check,build_mode)
% opt = geom.opt;
% if ~logical(getOptField(opt,'cmap',false))
%     error('buildLaplacePeanutBigSparse:UnsupportedOption', ...
%         'opt.use_big_sparse=1 requires opt.cmap=1.');
% end
% if getOptField(opt,'N_peanut',0) <= 0
%     error('buildLaplacePeanutBigSparse:UnsupportedOption', ...
%         'opt.use_big_sparse=1 requires opt.N_peanut > 0.');
% end
% if ~isequal(geom.rcheck,geom.rvec_out)
%     error('buildLaplacePeanutBigSparse:UnsupportedGrid', ...
%         'The Laplace big-sparse matvec is only supported on the solve grid.');
% end
% if N_check ~= round(N_check)
%     error('buildLaplacePeanutBigSparse:BadGridSize', ...
%         'rcheck length must be divisible by the number of particles.');
% end
% if strcmp(build_mode,'streaming') && ...
%         logical(getOptField(opt,'get_bndry_field',false))
%     error('buildLaplacePeanutBigSparse:StreamingBoundaryUnsupported', ...
%         ['lap_big_sparse_build_mode=''streaming'' is only supported when ', ...
%          'opt.get_bndry_field=0. Use ''auto'' or ''precomputed'' for ', ...
%          'boundary postprocessing.']);
% end
% if ~isfield(geom,'rout_base_f') || isempty(geom.rout_base_f)
%     error('buildLaplacePeanutBigSparse:MissingFineOutputGrid', ...
%         'geom.rout_base_f is required for Laplace big-sparse builds.');
% end
% end
% 
% function stats = initLaplaceBigSparseStats(n_pairs,N_c,N_check, ...
%     ram_estimate,build_mode)
% counts = ram_estimate.counts;
% stats = struct();
% stats.requested = false;
% stats.active = false;
% stats.backend = 'global_block_sparse';
% stats.reason = 'not_requested';
% stats.build_mode = build_mode;
% stats.n_pairs = n_pairs;
% stats.N_c = N_c;
% stats.N_check = N_check;
% stats.used_pair_cache = false;
% stats.rotations_used = false;
% stats.source_correction = 'factored_structured';
% stats.field_correction = 'factored_sparse';
% stats.local_pair_nonp_entries = counts.pair_nonp/max(1,n_pairs);
% stats.local_u_entries = counts.u/max(1,n_pairs);
% stats.local_qv_entries = counts.qv/max(1,n_pairs);
% stats.nnz_u = 0;
% stats.nnz_u_cross = 0;
% stats.nnz_u_peanut = 0;
% stats.nnz_u_qv = 0;
% stats.nnz_qv = 0;
% stats.nnz_qv_canon = 0;
% stats.nnz_pair_nonp = 0;
% stats.nnz_pair_nonp_canon = 0;
% stats.nnz_source_scatter = 0;
% stats.nnz_qv_scatter = 0;
% stats.estimated_sparse_bytes = ram_estimate.estimated_sparse_bytes;
% stats.estimated_auxiliary_bytes = ram_estimate.estimated_auxiliary_bytes;
% stats.estimated_build_bytes = ram_estimate.estimated_build_bytes;
% stats.estimated_peak_bytes = ram_estimate.estimated_peak_bytes;
% stats.estimated_sparse_MB = ram_estimate.estimated_sparse_MB;
% stats.estimated_auxiliary_MB = ram_estimate.estimated_auxiliary_MB;
% stats.estimated_build_MB = ram_estimate.estimated_build_MB;
% stats.estimated_peak_MB = ram_estimate.estimated_peak_MB;
% stats.max_build_bytes = inf;
% stats.max_build_MB = inf;
% stats.max_build_ram_fraction = [];
% stats.available_ram_bytes = nan;
% stats.available_ram_MB = nan;
% stats.guardrail_scope = 'estimated_peak_bytes';
% stats.guardrail_checked_bytes = ram_estimate.estimated_peak_bytes;
% stats.guardrail_checked_MB = ram_estimate.estimated_peak_MB;
% stats.build_time = 0;
% end
% 
% function build_limit = resolveLaplaceBigSparseBuildMemoryLimit(opt)
% explicit_bytes = getOptField(opt,'lap_big_sparse_max_build_bytes', ...
%     getOptField(opt,'big_sparse_max_build_bytes',inf));
% ram_fraction = getOptField(opt,'lap_big_sparse_max_build_ram_fraction',[]);
% if isempty(explicit_bytes)
%     explicit_bytes = inf;
% end
% 
% if isempty(ram_fraction) || (isnumeric(ram_fraction) && isnan(ram_fraction))
%     available_bytes = nan;
%     max_build_bytes = explicit_bytes;
% else
%     validateattributes(ram_fraction,{'numeric'}, ...
%         {'scalar','real','positive','<=',1},mfilename, ...
%         'opt.lap_big_sparse_max_build_ram_fraction');
%     available_bytes = getLinuxMemAvailableBytes();
%     max_build_bytes = min(explicit_bytes,ram_fraction*available_bytes);
% end
% 
% build_limit = struct('max_build_bytes',max_build_bytes, ...
%     'ram_fraction',ram_fraction, ...
%     'available_ram_bytes',available_bytes, ...
%     'scope','estimated_peak_bytes');
% end
% 
% function bytes = getLinuxMemAvailableBytes()
% meminfo_file = '/proc/meminfo';
% if exist(meminfo_file,'file') ~= 2
%     error('buildLaplacePeanutBigSparse:MemAvailableUnsupported', ...
%         ['opt.lap_big_sparse_max_build_ram_fraction requires ', ...
%          '/proc/meminfo. Set opt.lap_big_sparse_max_build_bytes directly ', ...
%          'on this platform. The guardrail is checked against ', ...
%          'estimated_peak_bytes.']);
% end
% 
% txt = fileread(meminfo_file);
% tok = regexp(txt,'MemAvailable:\s+(\d+)\s+kB','tokens','once');
% if isempty(tok)
%     error('buildLaplacePeanutBigSparse:MemAvailableMissing', ...
%         'Could not read MemAvailable from /proc/meminfo.');
% end
% bytes = str2double(tok{1})*1024;
% end
% 
% function big_sparse = emptyBigSparse(pair_total_rows,n_coarse,n_u_rows, ...
%     n_qv_rows,ram_estimate,P_pair)
% big_sparse = struct();
% big_sparse.matrix_plan = ram_estimate.matrix_plan;
% big_sparse.ram_estimate = ram_estimate;
% big_sparse.P_pair = P_pair;
% big_sparse.source_scatter_rows = zeros(0,1);
% big_sparse.qv_scatter_rows = zeros(0,1);
% big_sparse.pair_map_mode = 'direct';
% big_sparse.M_pair_nonp = sparse(pair_total_rows,n_coarse);
% big_sparse.M_u_cross = sparse(n_u_rows,n_coarse);
% big_sparse.M_u_peanut = sparse(n_u_rows,pair_total_rows);
% big_sparse.M_u_qv = sparse(n_u_rows,n_qv_rows);
% big_sparse.M_pair_qv = sparse(n_qv_rows,n_coarse);
% end
% 
% function pair_cache = emptyPairCacheForMode(build_mode,opt)
% pair_cache = initLaplacePairCache();
% pair_cache.enabled = logical(getOptField(opt,'reuse_pair_basis_by_sep',false));
% pair_cache.stats.branch = ternary(pair_cache.enabled,'canonical_group','per_pair');
% pair_cache.stats.n_pairs = 0;
% pair_cache.stats.n_groups = 0;
% pair_cache.stats.used_parallel = false;
% pair_cache.stats.build_mode = build_mode;
% end
% 
% function entries = initSparseBuilders(P,N_c,N_check,n_pairs,opt)
% chunk_pairs = max(1,round(getOptField(opt,'lap_big_sparse_chunk_pairs',8)));
% chunk_pairs = min(max(1,n_pairs),chunk_pairs);
% n_coarse = P*N_c;
% pair_rows = 2*N_c;
% use_factored_pair_maps = logical(getOptField(opt, ...
%     'reuse_pair_basis_by_sep',false));
% 
% entries = struct();
% entries.pair_map_mode = ternary(use_factored_pair_maps,'rotated_factored', ...
%     'direct');
% if use_factored_pair_maps
%     entries.pair_nonp = newSparseBuilder(0,n_pairs*pair_rows,n_coarse);
%     entries.pair_nonp_canon = newSparseBuilder( ...
%         chunk_pairs*pair_rows*pair_rows,n_pairs*pair_rows, ...
%         n_pairs*pair_rows);
% else
%     entries.pair_nonp = newSparseBuilder( ...
%         chunk_pairs*pair_rows*(2*N_c),n_pairs*pair_rows,n_coarse);
%     entries.pair_nonp_canon = newSparseBuilder(0,n_pairs*pair_rows, ...
%         n_pairs*pair_rows);
% end
% entries.u_cross = newSparseBuilder( ...
%     chunk_pairs*(2*N_check)*(2*N_c),P*N_check,n_coarse);
% entries.u_peanut = newSparseBuilder( ...
%     chunk_pairs*(2*N_check)*pair_rows,P*N_check,n_pairs*pair_rows);
% entries.u_qv = newSparseBuilder( ...
%     chunk_pairs*(2*N_check)*2,P*N_check,2*n_pairs);
% if use_factored_pair_maps
%     entries.qv = newSparseBuilder(0,2*n_pairs,n_coarse);
%     entries.qv_canon = newSparseBuilder( ...
%         chunk_pairs*2*pair_rows,2*n_pairs,n_pairs*pair_rows);
% else
%     entries.qv = newSparseBuilder( ...
%         chunk_pairs*2*(2*N_c),2*n_pairs,n_coarse);
%     entries.qv_canon = newSparseBuilder(0,2*n_pairs,n_pairs*pair_rows);
% end
% end
% 
% function block = newSparseBuilder(capacity,n_rows,n_cols)
% capacity = max(0,round(capacity));
% block = struct();
% block.rows = zeros(capacity,1);
% block.cols = zeros(capacity,1);
% block.vals = zeros(capacity,1);
% block.next = 1;
% block.S = sparse(n_rows,n_cols);
% end
% 
% function [entries,pair_cache,stats] = appendPrecomputedBlocks(entries, ...
%     geom,basis,P_pair,stats)
% pairs = geom.pairs;
% q = geom.q(:);
% opt = geom.opt;
% N_c = opt.N_c;
% n_pairs = size(pairs,1);
% N_check = numel(geom.rcheck)/numel(q);
% rout_base_c = geom.rcheck(1:N_check)-q(1);
% project_charge = logical(getOptField(opt,'project_charge',false));
% 
% if ~isfield(basis,'pair_cache') || isempty(basis.pair_cache)
%     error('buildLaplacePeanutBigSparse:MissingPairCache', ...
%         'precomputed mode requires basis.pair_cache.');
% end
% pair_cache = basis.pair_cache;
% use_pair_cache = isfield(pair_cache,'enabled') && pair_cache.enabled;
% stats.used_pair_cache = use_pair_cache;
% stats.rotations_used = use_pair_cache;
% 
% for row = 1:n_pairs
%     i = pairs(row,1);
%     j = pairs(row,2);
%     meta = pair_cache.meta(row);
%     if use_pair_cache
%         group = pair_cache.groups(meta.group_id);
%         [Ucross,Ec,Lr] = getActualCollocFactors(meta,q,geom.rbase_in_c, ...
%             rout_base_c);
%         entries = appendFactoredCanonicalPairBlocks(entries,pairs,row, ...
%             N_c,N_check,numel(q),group,Ucross,Ec,Lr,P_pair, ...
%             project_charge);
%     else
%         C_nonp = basis.Cmap{i,j};
%         QV = basis.Cmap_QV{i,j};
%         [Ucross,Ec,Lr] = getActualCollocFactors(meta,q,geom.rbase_in_c, ...
%             rout_base_c);
%         entries = appendPairBlocks(entries,pairs,row,N_c,N_check,numel(q), ...
%             C_nonp,QV,Ucross,Ec,Lr,P_pair,project_charge);
%     end
% end
% pair_cache.stats.build_mode = stats.build_mode;
% end
% 
% function [entries,pair_cache,stats] = appendStreamingBlocks(entries, ...
%     geom,P_pair,stats)
% q = geom.q(:);
% pairs = geom.pairs;
% opt = geom.opt;
% N_c = opt.N_c;
% N_f = opt.N_f;
% n_pairs = size(pairs,1);
% N_check = numel(geom.rcheck)/numel(q);
% rout_base_c = geom.rcheck(1:N_check)-q(1);
% project_charge = logical(getOptField(opt,'project_charge',false));
% reuse_pair_basis = logical(getOptField(opt,'reuse_pair_basis_by_sep',false));
% show_counter = logical(getOptField(opt,'show_counter',false));
% 
% pair_cache = initLaplacePairCache();
% pair_cache.enabled = reuse_pair_basis;
% pair_cache.check_rotations = false;
% pair_cache.shared_sep_tol = getOptField(opt,'shared_sep_tol',0);
% pair_cache.rout_base_f = geom.rout_base_f(:);
% pair_cache.meta = buildLaplacePairMeta(q,pairs,N_c,N_f,opt);
% pair_cache.stats.n_pairs = n_pairs;
% pair_cache.stats.used_parallel = false;
% pair_cache.stats.build_mode = stats.build_mode;
% stats.used_pair_cache = reuse_pair_basis;
% stats.rotations_used = reuse_pair_basis;
% 
% if ~reuse_pair_basis
%     pair_cache.stats.branch = 'per_pair';
%     for row = 1:n_pairs
%         group = buildLaplacePairGroup([],row,pair_cache.meta(row).sep, ...
%             q,geom.rbase_in_c,geom.rbase_in_f,geom.rimage_vec, ...
%             geom.refine,pairs,opt,geom.rout_base_f,project_charge,false);
%         [Ucross,Ec,Lr] = buildLaplaceActualPairCollocFactors( ...
%             pair_cache.meta(row),q,geom.rbase_in_c,rout_base_c);
%         entries = appendPairBlocks(entries,pairs,row,N_c,N_check,numel(q), ...
%             group.Cmap,group.Cmap_QV,Ucross,Ec,Lr,P_pair, ...
%             project_charge);
%         if show_counter
%             fprintf('buildLaplacePeanutBigSparse: streamed pair %d/%d\n', ...
%                 row,n_pairs);
%         end
%     end
%     return
% end
% 
% [group_id,group_sep,rep_rows] = groupLaplacePairSeparations( ...
%     pair_cache.meta,pair_cache.shared_sep_tol);
% pair_cache.group_id = group_id;
% pair_cache.group_sep = group_sep;
% pair_cache.representative_rows = rep_rows;
% pair_cache.n_groups = numel(group_sep);
% pair_cache.stats.branch = 'canonical_group';
% pair_cache.stats.n_groups = pair_cache.n_groups;
% 
% for row = 1:n_pairs
%     gid = group_id(row);
%     pair_cache.meta(row).group_id = gid;
%     pair_cache.meta(row).sep = group_sep(gid);
% end
% 
% for gg = 1:pair_cache.n_groups
%     group = buildLaplacePairGroup(gg,rep_rows(gg),group_sep(gg), ...
%         q,geom.rbase_in_c,geom.rbase_in_f,geom.rimage_vec, ...
%         geom.refine,pairs,opt,geom.rout_base_f,project_charge,true);
%     rows = find(group_id == gg).';
%     block_cache = containers.Map('KeyType','char','ValueType','any');
%     for row = rows
%         meta = pair_cache.meta(row);
%         cache_key = pairBlockCacheKey(meta.rot);
%         if isKey(block_cache,cache_key)
%             block = block_cache(cache_key);
%             Ucross = block.Ucross;
%             Ec = block.Ec;
%             Lr = block.Lr;
%         else
%             [Ucross,Ec,Lr] = buildLaplaceActualPairCollocFactors( ...
%                 meta,q,geom.rbase_in_c,rout_base_c);
%             block_cache(cache_key) = struct('Ucross',Ucross,'Ec',Ec, ...
%                 'Lr',Lr);
%         end
%         entries = appendFactoredCanonicalPairBlocks(entries,pairs,row,N_c, ...
%             N_check,numel(q),group,Ucross,Ec,Lr,P_pair, ...
%             project_charge);
%     end
%     if show_counter
%         fprintf(['buildLaplacePeanutBigSparse: streamed canonical group ', ...
%             '%d/%d covering %d pairs\n'],gg,pair_cache.n_groups,numel(rows));
%     end
% end
% end
% 
% function key = pairBlockCacheKey(rot)
% key = sprintf('%.16g_%.16g',real(rot),imag(rot));
% end
% 
% function [Ucross,Ec,Lr] = getActualCollocFactors(meta,q,rbase_in_c, ...
%     rout_base_c)
% if isfield(meta,'Ucross_colloc_actual') && ...
%         ~isempty(meta.Ucross_colloc_actual) && ...
%         isfield(meta,'Ec_colloc_actual') && ~isempty(meta.Ec_colloc_actual)
%     Ucross = meta.Ucross_colloc_actual;
%     Ec = meta.Ec_colloc_actual;
%     if isfield(meta,'Lr_colloc_actual') && ~isempty(meta.Lr_colloc_actual)
%         Lr = meta.Lr_colloc_actual;
%     else
%         Lr = buildLr(numel(rout_base_c));
%     end
% else
%     [Ucross,Ec,Lr] = buildLaplaceActualPairCollocFactors( ...
%         meta,q,rbase_in_c,rout_base_c);
% end
% end
% 
% function Lr = buildLr(n)
% Lr = [ones(n,1) zeros(n,1); zeros(n,1) ones(n,1)];
% end
% 
% function entries = appendPairBlocks(entries,pairs,row,N_c,N_check,P, ...
%     C_nonp,QV,Ucross,Ec,Lr,P_pair,project_charge)
% i = pairs(row,1);
% j = pairs(row,2);
% in_idx = pairCoarseInputIndices(i,j,N_c);
% global_in_idx = pairCoarseInputIndicesGlobal(i,j,N_c,P);
% pair_idx = (row-1)*(2*N_c)+1:row*(2*N_c);
% u_idx = pairOutputIndices(i,j,N_check);
% qv_idx = (row-1)*2+1:row*2;
% 
% entries.pair_nonp = appendDenseBlock(entries.pair_nonp,pair_idx, ...
%     global_in_idx,C_nonp);
% entries.u_cross = appendDenseBlock(entries.u_cross,u_idx,global_in_idx, ...
%     Ucross);
% entries.u_peanut = appendDenseBlock(entries.u_peanut,u_idx,pair_idx,Ec);
% if project_charge
%     entries.u_qv = appendDenseBlock(entries.u_qv,u_idx,qv_idx,Lr);
% end
% entries.qv = appendDenseBlock(entries.qv,qv_idx,global_in_idx,QV);
% 
% % in_idx is intentionally retained as a local consistency check.
% if numel(in_idx) ~= 2*N_c
%     error('buildLaplacePeanutBigSparse:BadLocalIndex', ...
%         'Internal pair input index size mismatch.');
% end
% end
% 
% function entries = appendFactoredCanonicalPairBlocks(entries,pairs,row,N_c, ...
%     N_check,P,group,Ucross,Ec,Lr,P_pair,project_charge)
% i = pairs(row,1);
% j = pairs(row,2);
% global_in_idx = pairCoarseInputIndicesGlobal(i,j,N_c,P);
% pair_idx = (row-1)*(2*N_c)+1:row*(2*N_c);
% u_idx = pairOutputIndices(i,j,N_check);
% qv_idx = (row-1)*2+1:row*2;
% 
% entries.pair_nonp_canon = appendDenseBlock(entries.pair_nonp_canon, ...
%     pair_idx,pair_idx,group.Cmap);
% 
% entries.u_cross = appendDenseBlock(entries.u_cross,u_idx,global_in_idx, ...
%     Ucross);
% entries.u_peanut = appendDenseBlock(entries.u_peanut,u_idx,pair_idx,Ec);
% if project_charge
%     entries.u_qv = appendDenseBlock(entries.u_qv,u_idx,qv_idx,Lr);
% end
% entries.qv_canon = appendDenseBlock(entries.qv_canon,qv_idx,pair_idx, ...
%     group.Cmap_QV);
% 
% if size(P_pair,1) ~= 2*N_c
%     error('buildLaplacePeanutBigSparse:BadPairProjector', ...
%         'Internal pair projector size mismatch.');
% end
% end
% 
% function block = appendDenseBlock(block,row_idx,col_idx,A)
% row_idx = row_idx(:);
% col_idx = col_idx(:);
% nr = numel(row_idx);
% nc = numel(col_idx);
% n = nr*nc;
% if n == 0
%     return
% end
% if block.next+n-1 > numel(block.vals)
%     block = flushSparseBuilder(block);
% end
% if n > numel(block.vals)
%     block.S = block.S + sparse(repmat(row_idx,nc,1), ...
%         repelem(col_idx,nr),real(A(:)),size(block.S,1),size(block.S,2));
%     return
% end
% loc = block.next:block.next+n-1;
% block.rows(loc) = repmat(row_idx,nc,1);
% block.cols(loc) = repelem(col_idx,nr);
% block.vals(loc) = real(A(:));
% block.next = block.next + n;
% end
% 
% function block = flushSparseBuilder(block)
% n = block.next - 1;
% if n > 0
%     block.S = block.S + sparse(block.rows(1:n),block.cols(1:n), ...
%         block.vals(1:n),size(block.S,1),size(block.S,2));
%     block.next = 1;
% end
% end
% 
% function entries = flushAllSparseBuilders(entries)
% entries.pair_nonp = flushSparseBuilder(entries.pair_nonp);
% entries.pair_nonp_canon = flushSparseBuilder(entries.pair_nonp_canon);
% entries.u_cross = flushSparseBuilder(entries.u_cross);
% entries.u_peanut = flushSparseBuilder(entries.u_peanut);
% entries.u_qv = flushSparseBuilder(entries.u_qv);
% entries.qv = flushSparseBuilder(entries.qv);
% entries.qv_canon = flushSparseBuilder(entries.qv_canon);
% end
% 
% function P_pair = buildChargeProjectionPairMatrix(N_c,project_charge)
% if project_charge
%     P_body = eye(N_c) - ones(N_c,N_c)/N_c;
% else
%     P_body = eye(N_c);
% end
% Z = zeros(N_c,N_c);
% P_pair = [P_body Z; Z P_body];
% end
% 
% function rows = buildSourceScatterRows(pairs,N_c,P)
% n_pairs = size(pairs,1);
% rows = zeros(n_pairs*4*N_c,1);
% next = 1;
% n_coarse = P*N_c;
% for row = 1:n_pairs
%     i = pairs(row,1);
%     j = pairs(row,2);
%     idx_i = (i-1)*N_c+1:i*N_c;
%     idx_j = (j-1)*N_c+1:j*N_c;
%     rows(next:next+4*N_c-1) = [idx_i idx_j ...
%         n_coarse+idx_i n_coarse+idx_j].';
%     next = next + 4*N_c;
% end
% end
% 
% function rows = buildQvScatterRows(pairs)
% n_pairs = size(pairs,1);
% rows = zeros(2*n_pairs,1);
% for row = 1:n_pairs
%     rows((row-1)*2+1:row*2) = pairs(row,:).';
% end
% end
% 
% function idx = pairCoarseInputIndicesGlobal(i,j,N_c,P)
% idx_i = (i-1)*N_c+1:i*N_c;
% idx_j = (j-1)*N_c+1:j*N_c;
% idx = [idx_i idx_j].';
% if max(idx) > P*N_c
%     error('buildLaplacePeanutBigSparse:BadGlobalIndex', ...
%         'Internal global input index out of range.');
% end
% end
% 
% function idx = pairCoarseInputIndices(i,j,N_c) %#ok<INUSD>
% idx = (1:2*N_c).';
% end
% 
% function idx = pairOutputIndices(i,j,N_check)
% idx_i = (i-1)*N_check+1:i*N_check;
% idx_j = (j-1)*N_check+1:j*N_check;
% idx = [idx_i idx_j].';
% end
% 
% function out = ternary(cond,a,b)
% if cond
%     out = a;
% else
%     out = b;
% end
% end
% 
% function test_buildLaplacePeanutBigSparse
% fprintf('buildLaplacePeanutBigSparse self-test: shared-pair spy plots\n');
% rng(220422);
% 
% [geom,basis] = buildSelfTestLaplaceData();
% tau = randn(numel(geom.rvec_out),1);
% res_ref = matvec_lap_peanut_enhanced(tau,geom,basis);
% 
% geom_pre = geom;
% geom_pre.opt.use_big_sparse = true;
% geom_pre.opt.lap_big_sparse_build_mode = 'precomputed';
% [big_pre,stats_pre,pair_cache_pre] = buildLaplacePeanutBigSparse( ...
%     geom_pre,basis);
% basis_pre = basis;
% basis_pre.big_sparse = big_pre;
% basis_pre.pair_cache = pair_cache_pre;
% geom_pre.pair_cache = pair_cache_pre;
% res_pre = matvec_lap_peanut_big_sparse(tau,geom_pre,basis_pre);
% 
% geom_stream = geom_pre;
% geom_stream.opt.lap_big_sparse_build_mode = 'streaming';
% basis_stream = basis;
% basis_stream.Upf = [];
% basis_stream.Ypf = [];
% basis_stream.DC_all = [];
% basis_stream.YC_all = [];
% basis_stream.Cmap = [];
% basis_stream.Cmap_QV = [];
% basis_stream.pair_cache = initLaplacePairCache();
% [big_stream,stats_stream,pair_cache_stream] = buildLaplacePeanutBigSparse( ...
%     geom_stream,basis_stream);
% basis_stream.big_sparse = big_stream;
% basis_stream.pair_cache = pair_cache_stream;
% geom_stream.pair_cache = pair_cache_stream;
% res_stream = matvec_lap_peanut_big_sparse(tau,geom_stream,basis_stream);
% 
% err_pre = relerrSelfTest(res_pre,res_ref);
% err_stream = relerrSelfTest(res_stream,res_ref);
% 
% fprintf('  P=%d, close pairs=%d, shared groups=%d\n', ...
%     numel(geom.q),stats_pre.n_pairs,pair_cache_pre.n_groups);
% fprintf(['  precomputed: build %.3fs, sparse %.2f MB, ', ...
%     'canonical pair nnz=%d, relerr %.3e\n'], ...
%     stats_pre.build_time,stats_pre.estimated_sparse_MB, ...
%     stats_pre.nnz_pair_nonp_canon,err_pre);
% fprintf(['  streaming:   build %.3fs, sparse %.2f MB, ', ...
%     'canonical pair nnz=%d, relerr %.3e\n'], ...
%     stats_stream.build_time,stats_stream.estimated_sparse_MB, ...
%     stats_stream.nnz_pair_nonp_canon,err_stream);
% 
% figure('Name','buildLaplacePeanutBigSparse spy self-test','Color','w');
% tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
% 
% nexttile;
% spy(big_pre.M_pair_nonp_canon);
% title('precomputed: M_pair_nonp_canon','Interpreter','none');
% 
% nexttile;
% spy(big_pre.M_u_cross);
% title('precomputed: M_u_cross','Interpreter','none');
% 
% nexttile;
% spy(big_pre.M_u_peanut);
% title('precomputed: M_u_peanut','Interpreter','none');
% 
% nexttile;
% spy(big_stream.M_pair_nonp_canon);
% title('streaming: M_pair_nonp_canon','Interpreter','none');
% 
% nexttile;
% spy(big_stream.M_u_cross);
% title('streaming: M_u_cross','Interpreter','none');
% 
% nexttile;
% spy(big_stream.M_u_peanut);
% title('streaming: M_u_peanut','Interpreter','none');
% 
%     function [geom,basis] = buildSelfTestLaplaceData()
%     R = 2;
%     gap = 2e-2;
% 
%     P = 50;
%     theta = 1.3;
%     q = exp(1i*theta)*grow_cluster(P,gap,2,2);
%     N_c = 20;
%     N_f = 28;
% 
%     opt = getLaplace2Dparams(P,R,N_c,N_f);
%     opt.N_peanut = 64;
%     opt.delta_pair = 0.2*R;
%     opt.cmap = 1;
%     opt.compress_cmap = 0;
%     opt.reuse_pair_basis_by_sep = true;
%     opt.shared_sep_tol = 1e-8;
%     opt.get_bndry_field = 0;
%     opt.use_fmm = false;
%     opt.show_counter = 0;
%     opt.visualise_sol = 0;
%     opt.visualise_grid = 0;
%     opt.gmres_verbose = 0;
% 
%     [rbase_in_c,rbase_out_c,rbase_in_f,rout_base_f,nout] = ...
%         buildSelfTestCircleGrids(R,opt);
% 
%     rvec_in_c = zeros(P*N_c,1);
%     rout = zeros(P*nout,1);
%     for k = 1:P
%         rvec_in_c((k-1)*N_c+1:k*N_c) = q(k)+rbase_in_c;
%         rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
%     end
% 
%     [~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
%     [UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_QV,pair_cache] = ...
%         getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f, ...
%         rbase_out_c,rimage_vec,refine,pairs,opt);
%     [U,Y] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout]);
% 
%     geom = struct();
%     geom.rbase_in_c = rbase_in_c;
%     geom.rbase_in_f = rbase_in_f;
%     geom.rout_base_f = rout_base_f;
%     geom.refine = refine;
%     geom.opt = opt;
%     geom.rvec_out = rout;
%     geom.rcheck = rout;
%     geom.q = q;
%     geom.pairs = pairs;
%     geom.rimage_vec = rimage_vec;
%     geom.rvec_in = rvec_in_c;
%     geom.pair_cache = pair_cache;
% 
%     basis = struct();
%     basis.U = U;
%     basis.Y = Y;
%     basis.Upf = UB_all;
%     basis.Ypf = YB_all;
%     basis.DC_all = UC_all;
%     basis.YC_all = YC_all;
%     basis.Cmap = Cmap;
%     basis.Cmap_QV = Cmap_QV;
%     basis.pair_cache = pair_cache;
%     basis.Nii = lapSLPmat(rbase_in_c,rbase_out_c);
%     end
% 
%     function [rbase_in_c,rbase_out_c,rbase_in_f,rout_base_f,nout] = ...
%             buildSelfTestCircleGrids(R,opt)
%     N_c = opt.N_c;
%     N_f = opt.N_f;
%     nout = ceil(opt.a_c*N_c);
% 
%     t = linspace(0,2*pi,N_c+1)';
%     t = t(1:end-1);
%     rbase_in_c = opt.Rp_c*(cos(t)+1i*sin(t));
% 
%     t = linspace(0,2*pi,nout+1)';
%     t = t(1:end-1);
%     rbase_out_c = R*(cos(t)+1i*sin(t));
% 
%     t = linspace(0,2*pi,N_f+1)';
%     t = t(1:end-1);
%     rbase_in_f = opt.Rp_f*(cos(t)+1i*sin(t));
% 
%     t = linspace(0,2*pi,ceil(opt.a_f*N_f)+1)';
%     t = t(1:end-1);
%     rout_base_f = R*(cos(t)+1i*sin(t));
%     end
% 
%     function e = relerrSelfTest(a,b)
%     e = norm(a-b,inf)/max(1,norm(b,inf));
%     end
% end
