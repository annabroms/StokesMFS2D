function result = test_getPairBasisStokes_parallel
%TEST_GETPAIRBASISSTOKES_PARALLEL Small parity checks for chunked Stokes pair precompute.

fprintf('--- test_getPairBasisStokes_parallel ---\n');

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

ensure_pool(2);

P = 10;
q = grow_cluster(P,1e-3,2);
base_opt = make_small_options(P);

result = struct();
result.cases = repmat(init_case_result(),0,1);
case_id = 0;
for get_bndry_field = [false true]
    case_id = case_id + 1;
    result.cases(case_id) = compare_case(q,base_opt,get_bndry_field);
end

assert(max([result.cases.max_pair_map_rel_diff]) < 5e-10, ...
    'Chunked parallel pair precompute changed Cmap or Cmap_FU.');
assert(max([result.cases.max_full_payload_rel_diff]) < 1e-4, ...
    'Chunked parallel pair precompute changed the stable full-payload operator.');
assert(max([result.cases.max_dense_actual_rel_diff]) < 5e-10, ...
    'Chunked parallel pair precompute changed actual dense collocation fields.');

fprintf('  PASSED\n');
end

function case_result = compare_case(q,base_opt,get_bndry_field)
[rbase_in_c,rbase_in_f,rbase_out_c,rimage_vec,refine,pairs,Lc,svd_opts,opt_pair] = ...
    build_pair_setup(q,base_opt,get_bndry_field);

opt_serial = opt_pair;
opt_serial.parallel_precomp = false;

opt_parallel = opt_pair;
opt_parallel.parallel_precomp = true;
opt_parallel.parallel_precomp_chunk_pairs = 2;

[Uf_s,Yf_s,Up_s,Yp_s,Cmap_s,Cmap_FU_s,cache_s] = getPairBasisStokes( ...
    q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt_serial,Lc,rbase_out_c,svd_opts);
[Uf_p,Yf_p,Up_p,Yp_p,Cmap_p,Cmap_FU_p,cache_p] = getPairBasisStokes( ...
    q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt_parallel,Lc,rbase_out_c,svd_opts);

case_result = init_case_result();
case_result.get_bndry_field = logical(get_bndry_field);
case_result.payload_mode = cache_p.stats.payload_mode;
case_result.n_pairs = cache_p.stats.n_pairs;
case_result.parallel_backend = cache_p.stats.parallel_backend;
case_result.chunk_pairs = cache_p.stats.chunk_pairs;
case_result.n_tasks = cache_p.stats.n_tasks;
case_result.pool_size = cache_p.stats.pool_size;
case_result.max_pair_map_rel_diff = max( ...
    compare_cell_grid(Cmap_p,Cmap_s),compare_cell_grid(Cmap_FU_p,Cmap_FU_s));
case_result.max_full_payload_rel_diff = compare_full_payload_operator( ...
    Uf_p,Yf_p,Uf_s,Yf_s);
case_result.max_dense_actual_rel_diff = compare_dense_meta(cache_p,cache_s);

assert(strcmp(cache_p.stats.parallel_backend,'parfor_chunked'), ...
    'Parallel pair precompute did not use the chunked backend.');
assert(cache_p.stats.n_tasks >= 1,'Parallel pair precompute recorded no tasks.');
assert(cache_p.stats.chunk_pairs >= 1,'Parallel pair precompute recorded no chunk size.');

fprintf(['  get_bndry_field=%d payload=%s pairs=%d chunk=%d tasks=%d ', ...
    'pair_map_diff=%.3e full_diff=%.3e dense_diff=%.3e\n'], ...
    get_bndry_field,case_result.payload_mode,case_result.n_pairs, ...
    case_result.chunk_pairs,case_result.n_tasks, ...
    case_result.max_pair_map_rel_diff, ...
    case_result.max_full_payload_rel_diff, ...
    case_result.max_dense_actual_rel_diff);
end

function case_result = init_case_result()
case_result = struct( ...
    'get_bndry_field',false, ...
    'payload_mode','', ...
    'n_pairs',0, ...
    'parallel_backend','', ...
    'chunk_pairs',0, ...
    'n_tasks',0, ...
    'pool_size',0, ...
    'max_pair_map_rel_diff',0, ...
    'max_full_payload_rel_diff',0, ...
    'max_dense_actual_rel_diff',0);
end

function [rbase_in_c,rbase_in_f,rbase_out_c,rimage_vec,refine,pairs,Lc,svd_opts,opt_pair] = ...
    build_pair_setup(q,base_opt,get_bndry_field)
q = q(:);
P = numel(q);
rad = 1;
N_c = base_opt.N_c;
N_f = base_opt.N_f;
a_c = base_opt.a_c;
Rp_c = base_opt.Rp_c;
Rp_f = base_opt.Rp_f;

tout_c_all = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c_all(1:end-1)';
rbase_out_c = rad*(cos(tout_c)+1i*sin(tout_c));
tin_c = linspace(0,2*pi,N_c+1)';
tin_c = tin_c(1:end-1);
rbase_in_c = Rp_c*cos(tin_c) + 1i*Rp_c*sin(tin_c);
tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
rbase_in_f = Rp_f*cos(tin_f) + 1i*Rp_f*sin(tin_f);

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,base_opt);

svd_opts = struct('column_weight',false,'left_weight',false);
[~,~,Lc_cells] = getSelfPseudoMobilityStokes( ...
    1,q,rbase_in_c,rbase_out_c,[],[0,ceil(a_c*N_c)],svd_opts);
Lc = Lc_cells{1};

opt_pair = base_opt;
opt_pair.P = P;
opt_pair.get_bndry_field = logical(get_bndry_field);
opt_pair.project_force = true;
opt_pair.project = true;
opt_pair.pair_basis_debug = false;
opt_pair.rad = ones(P,1);
end

function opt = make_small_options(P)
opt = get2Dparams(P,12,18);
opt.delta_pair = 0.2;
opt.N_peanut = 48;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.parallel_precomp_chunk_pairs = 2;
opt.mob_big_sparse_chunk_pairs = 2;
opt.use_fmm = false;
end

function value = compare_cell_grid(a,b)
if isempty(a) && isempty(b)
    value = 0;
    return
end

assert(isequal(size(a),size(b)), ...
    'Cell output size mismatch in pair-basis comparison.');
value = 0;
for kk = 1:numel(a)
    value = max(value,relerr_optional(a{kk},b{kk}));
end
end

function value = compare_dense_meta(cache_a,cache_b)
value = 0;
fields = {'Upair_colloc_actual','Ucross_colloc_actual','Ecolloc_actual'};
assert(numel(cache_a.meta) == numel(cache_b.meta), ...
    'pair_cache.meta size mismatch in pair-basis comparison.');
for row = 1:numel(cache_a.meta)
    for ff = 1:numel(fields)
        value = max(value,relerr_optional( ...
            cache_a.meta(row).(fields{ff}),cache_b.meta(row).(fields{ff})));
    end
end
end

function value = compare_full_payload_operator(Uf_a,Yf_a,Uf_b,Yf_b)
if isempty(Uf_a) && isempty(Uf_b) && isempty(Yf_a) && isempty(Yf_b)
    value = 0;
    return
end

assert(isequal(size(Uf_a),size(Uf_b)) && isequal(size(Yf_a),size(Yf_b)), ...
    'Full payload cell output size mismatch in pair-basis comparison.');
value = 0;
for kk = 1:numel(Uf_a)
    if isempty(Uf_a{kk}) && isempty(Uf_b{kk}) && ...
            isempty(Yf_a{kk}) && isempty(Yf_b{kk})
        continue
    end
    assert(~isempty(Uf_a{kk}) && ~isempty(Uf_b{kk}) && ...
        ~isempty(Yf_a{kk}) && ~isempty(Yf_b{kk}), ...
        'Only one side of a full payload pair basis factor is empty.');
    value = max(value,relerr_optional(Yf_a{kk}*Uf_a{kk},Yf_b{kk}*Uf_b{kk}));
end
end

function value = relerr_optional(a,b)
if isempty(a) && isempty(b)
    value = 0;
    return
end

assert(~isempty(a) && ~isempty(b), ...
    'Only one side of a pair-basis output is empty.');
value = norm(a(:)-b(:),inf) / max(1,norm(b(:),inf));
end

function ensure_pool(pool_size)
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers ~= pool_size
    delete(pool);
    pool = [];
end
if isempty(pool)
    parpool('local',pool_size);
end
end
