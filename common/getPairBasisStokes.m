function [Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair_cache] = getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,rout_base_c,svd_opts)
%GETPAIRBASISSTOKES Build Stokes pair-basis pseudoinverse factors.
%
% Syntax:
%   [Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair_cache] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt)
%   [Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair_cache] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc)
%   [Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair_cache] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,rout_base_c)
%   [Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair_cache] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,rout_base_c,svd_opts)
%
% If opt.use_dense is enabled and rout_base_c is supplied, pair_cache.meta
% also stores the actual dense Stokeslet blocks on the solver collocation
% grid so the transform helpers can skip direct summation of Stokeslets.
%
% Anna Broms, Mar 2026

q = q(:);
P = opt.P;

N_peanut = opt.N_peanut; % number of peanut collocation nodes
use_pair_map = logical(opt.cmap); % use coarse-to-coarse compression?
get_bndry_field = logical(opt.get_bndry_field); % evaluate flow field on boundary? 
% If so, the mapping to recover fine sources must be stored
reuse_pair_basis = logical(opt.reuse_pair_basis_by_sep); % determine one compression per separation?
parallel_precomp = logical(opt.parallel_precomp); % parallelise pair precomputation?
shared_sep_tol = opt.shared_sep_tol; % tolerance for what is considered the same separation.
use_dense = logical(opt.use_dense); 

%Draw pair discretisation for debugging purpose?
debug = logical(opt.pair_basis_debug);
%Show progress for building pair basis?
show_counter = opt.show_counter; 
% Need fine sources explicitly?
self_correct = logical(getOptField(opt,'self_correct',false));
need_explicit_pair_sources = needs_explicit_pair_sources( ...
    N_peanut,use_pair_map,get_bndry_field,self_correct);
payload_mode = choose_pair_payload_mode( ...
    N_peanut,use_pair_map,get_bndry_field,self_correct,debug);
store_full_pair_payload = strcmp(payload_mode,'full');


pair_rad = getOptField(opt,'rad',1);
if numel(pair_rad) > 1
    pair_rad = pair_rad(1);
end

if nargin < 9
    rout_base_c = [];
end
if nargin < 10 || isempty(svd_opts)
    svd_opts = struct();
end
if ~isempty(rout_base_c)
    rout_base_c = rout_base_c(:);
end



if reuse_pair_basis
    Uf = [];
    Yf = [];
    Up = [];
    Yp = [];
    Cmap = [];
    Cmap_FU = [];
else
    [Uf,Yf,Up,Yp,Cmap,Cmap_FU] = init_outputs( ...
        P,N_peanut,use_pair_map,store_full_pair_payload);
end

pair_cache = init_pair_cache();
pair_cache.enabled = reuse_pair_basis;
pair_cache.shared_sep_tol = shared_sep_tol;
pair_cache.stats.requested_parallel = parallel_precomp;
pair_cache.stats.payload_mode = payload_mode;
pair_cache.stats.needs_explicit_pair_sources = need_explicit_pair_sources;

if isempty(pairs)
    return
end

N_f = opt.N_f;
pair_cache.meta = build_pair_meta(q,pairs,numel(rbase_in_c),N_f,opt);
pair_cache.stats.n_pairs = size(pairs,1);

if ~reuse_pair_basis
    total_pairs = size(pairs,1);
    pair_cache.stats.branch = 'per_pair';
    use_parallel_pairs = parallel_precomp && total_pairs > 1 && ~debug;
    [use_parallel_pairs,pool_size] = resolve_parallel_precomp(use_parallel_pairs, ...
        'getPairBasisStokes');
    pair_cache.stats.used_parallel = use_parallel_pairs;
    pair_cache.stats.pool_size = pool_size;

    if use_parallel_pairs
        if show_counter
            fprintf(['getPairBasisStokes: streamed parallel pair build for ', ...
                '%d pairs (%s payload)\n'], total_pairs,payload_mode);
        end
        max_inflight = max(1,pool_size);
        pair_cache.stats.max_inflight = max_inflight;
        pair_cache.stats.parallel_backend = 'parfeval_stream';
        [Uf,Yf,Up,Yp,Cmap,Cmap_FU,pool_size] = build_pairs_parallel_streamed( ...
            q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc, ...
            svd_opts,payload_mode,store_full_pair_payload,Uf,Yf,Up,Yp, ...
            Cmap,Cmap_FU,max_inflight,show_counter);
        pair_cache.stats.pool_size = pool_size;
        if show_counter
            fprintf('getPairBasisStokes: finished parallel pair build for %d pairs\n', ...
                total_pairs);
        end
    else
        pair_cache.stats.parallel_backend = 'serial';
        processed_pairs = 0;
        for ii = 1:total_pairs
            i = pairs(ii,1);
            p2 = pairs(ii,2);
            pair = build_pair_data(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs, ...
                opt,Lc,ii,debug,svd_opts,false,payload_mode);
            [Uf,Yf,Up,Yp,Cmap,Cmap_FU] = store_pair_payload( ...
                Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair,pairs,ii,N_peanut, ...
                use_pair_map,store_full_pair_payload);

            processed_pairs = processed_pairs + 1;
            if show_counter
                fprintf('getPairBasisStokes: processed pair %d/%d (%d,%d)\n', ...
                    processed_pairs,total_pairs,i,p2);
            end
        end
    end

    if use_dense && store_full_pair_payload && ~isempty(rout_base_c)
        pair_cache = populate_actual_dense_pair_fields(pair_cache,q,rbase_in_c,rbase_in_f, ...
            rout_base_c,rimage_pairs,pairs);
    end

    return
end

[group_id,group_sep,rep_rows] = group_pair_separations(pair_cache.meta,shared_sep_tol);
n_groups = numel(group_sep);
pair_cache.n_groups = n_groups;
pair_cache.group_id = group_id;
pair_cache.group_sep = group_sep;
pair_cache.representative_rows = rep_rows;
total_pairs = size(pairs,1);
pair_cache.stats.branch = 'canonical_group';
pair_cache.stats.n_groups = n_groups;
use_parallel_groups = parallel_precomp && n_groups > 1;
[use_parallel_groups,pool_size] = resolve_parallel_precomp(use_parallel_groups, ...
    'getPairBasisStokes');
pair_cache.stats.used_parallel = use_parallel_groups;
pair_cache.stats.pool_size = pool_size;

if use_parallel_groups
    pair_cache.stats.parallel_backend = 'parfor';
    pair_cache.stats.max_inflight = max(1,pool_size);
    if show_counter
        fprintf(['getPairBasisStokes: parallel canonical build for %d groups ', ...
            'covering %d pairs (%s payload)\n'], n_groups,total_pairs,payload_mode);
    end
    groups = repmat(init_pair_group(),n_groups,1);
    parfor gg = 1:n_groups
        groups(gg) = build_pair_group(gg,rep_rows(gg),group_sep(gg), ...
            q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc, ...
            store_full_pair_payload,false,svd_opts,payload_mode);
    end
    pair_cache.groups = groups;
    pair_cache.stats.pool_size = get_parallel_pool_size();
else
    pair_cache.stats.parallel_backend = 'serial';
    pair_cache.groups = repmat(init_pair_group(),n_groups,1);
    covered_pairs = 0;
    for gg = 1:n_groups
        pair_cache.groups(gg) = build_pair_group(gg,rep_rows(gg),group_sep(gg), ...
            q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc, ...
            store_full_pair_payload,false,svd_opts,payload_mode);

        if show_counter
            rep_i = pairs(rep_rows(gg),1);
            rep_j = pairs(rep_rows(gg),2);
            n_rep = sum(group_id == gg);
            covered_pairs = covered_pairs + n_rep;
            fprintf(['getPairBasisStokes: processed canonical group %d/%d ', ...
                'from pair (%d,%d), pp_sep = %.3g, covers %d pairs -> %d/%d pairs covered\n'], ...
                gg,n_groups,rep_i,rep_j,group_sep(gg)-2*pair_rad,n_rep,covered_pairs,total_pairs);
        end
    end
end

if debug
    % In reuse mode, canonical groups are only built for representative
    % pairs. Run the debug LS checks on every actual pair as requested.
    for row = 1:size(pairs,1)
        build_pair_data(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs, ...
            opt,Lc,row,debug,svd_opts,false,'full');
    end
end

for row = 1:size(pairs,1)
    gid = group_id(row);
    pair_cache.meta(row).group_id = gid;
    pair_cache.meta(row).sep = group_sep(gid);
end

if use_dense && store_full_pair_payload && ~isempty(rout_base_c)
    pair_cache = populate_actual_dense_pair_fields(pair_cache,q,rbase_in_c,rbase_in_f, ...
        rout_base_c,rimage_pairs,pairs);
end

if show_counter
    fprintf('getPairBasisStokes: reused %d canonical groups for %d pairs\n', ...
        n_groups,size(pairs,1));
end

end

function [Uf,Yf,Up,Yp,Cmap,Cmap_FU] = init_outputs(P,N_peanut,use_pair_map,need_explicit_pair_sources)
if need_explicit_pair_sources
    Uf = cell(P);
    Yf = cell(P);
else
    Uf = [];
    Yf = [];
end

if N_peanut && ~use_pair_map
    Up = cell(P);
    Yp = cell(P);
else
    Up = [];
    Yp = [];
end

if N_peanut && use_pair_map
    Cmap = cell(P);
else
    Cmap = [];
end

if use_pair_map
    Cmap_FU = cell(P);
else
    Cmap_FU = [];
end
end

function tf = needs_explicit_pair_sources(N_peanut,use_pair_map,get_bndry_field,self_correct)
tf = (N_peanut == 0) || ~use_pair_map || get_bndry_field || ~self_correct;
end

function payload_mode = choose_pair_payload_mode(N_peanut,use_pair_map, ...
    get_bndry_field,self_correct,debug)
if (N_peanut > 0) && use_pair_map && self_correct && ~get_bndry_field && ~debug
    payload_mode = 'maps_only';
else
    payload_mode = 'full';
end
end

function pair_cache = init_pair_cache()
pair_cache = struct();
pair_cache.enabled = false;
pair_cache.shared_sep_tol = [];
pair_cache.meta = repmat(struct('i',[],'j',[],'group_id',[],'sep',[], ...
    'mid',[],'rot',[],'phase_c',[],'phase_c_inv',[], ...
    'phase_f',[],'phase_f_inv',[], ...
    'Upair_colloc_actual',[],'Ucross_colloc_actual',[], ...
    'Ecolloc_actual',[]),0,1);
pair_cache.groups = repmat(init_pair_group(),0,1);
pair_cache.group_id = zeros(0,1);
pair_cache.group_sep = zeros(0,1);
pair_cache.representative_rows = zeros(0,1);
pair_cache.n_groups = 0;
pair_cache.stats = init_pair_stats();
end

function stats = init_pair_stats()
stats = struct('requested_parallel',false,'used_parallel',false, ...
    'branch','','n_pairs',0,'n_groups',0,'pool_size',0, ...
    'payload_mode','full','parallel_backend','none','max_inflight',0, ...
    'needs_explicit_pair_sources',false);
end

function [use_parallel,pool_size] = resolve_parallel_precomp(requested_parallel,caller_name)
use_parallel = requested_parallel;
pool_size = 0;
if ~requested_parallel
    return
end

if isempty(ver('parallel')) || ~license('test','Distrib_Computing_Toolbox') || ...
        exist('gcp','file') ~= 2
    error([caller_name ':ParallelToolboxRequired'], ...
        ['opt.parallel_precomp requires Parallel Computing Toolbox. ', ...
         'Open a pool before benchmarking warm-pool speedups if you want ', ...
         'to exclude startup overhead.']);
end

pool = gcp('nocreate');
if isempty(pool)
    pool = gcp();
end
pool_size = pool.NumWorkers;
end

function pool_size = get_parallel_pool_size()
pool_size = 0;
pool = gcp('nocreate');
if ~isempty(pool)
    pool_size = pool.NumWorkers;
end
end

function [Uf,Yf,Up,Yp,Cmap,Cmap_FU,pool_size] = build_pairs_parallel_streamed( ...
    q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,svd_opts, ...
    payload_mode,store_full_pair_payload,Uf,Yf,Up,Yp,Cmap,Cmap_FU, ...
    max_inflight,show_counter)

total_pairs = size(pairs,1);
ctx_data = struct();
ctx_data.q = q;
ctx_data.rbase_in_c = rbase_in_c;
ctx_data.rbase_in_f = rbase_in_f;
ctx_data.rimage_pairs = rimage_pairs;
ctx_data.refine = refine;
ctx_data.pairs = pairs;
ctx_data.opt = opt;
ctx_data.Lc = Lc;
ctx_data.svd_opts = svd_opts;
ctx = parallel.pool.Constant(ctx_data);

futures = parallel.FevalFuture.empty(0,1);
next_row = 1;
processed_pairs = 0;

[futures,next_row] = submit_pair_futures( ...
    futures,next_row,total_pairs,max_inflight,ctx,payload_mode);

while processed_pairs < total_pairs
    [completed_idx,row,pair] = fetchNext(futures);
    futures(completed_idx) = [];

    [Uf,Yf,Up,Yp,Cmap,Cmap_FU] = store_pair_payload( ...
        Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair,pairs,row,opt.N_peanut, ...
        logical(opt.cmap),store_full_pair_payload);

    processed_pairs = processed_pairs + 1;
    if show_counter
        i = pairs(row,1);
        p2 = pairs(row,2);
        fprintf('getPairBasisStokes: processed streamed pair %d/%d (%d,%d)\n', ...
            processed_pairs,total_pairs,i,p2);
    end

    [futures,next_row] = submit_pair_futures( ...
        futures,next_row,total_pairs,max_inflight,ctx,payload_mode);
end

clear futures ctx
pool_size = get_parallel_pool_size();
end

function [futures,next_row] = submit_pair_futures(futures,next_row, ...
    total_pairs,max_inflight,ctx,payload_mode)
while next_row <= total_pairs && numel(futures) < max_inflight
    futures(end+1) = parfeval(@build_pair_data_from_context, ...
        2,ctx,next_row,payload_mode);
    next_row = next_row + 1;
end
end

function [row,pair] = build_pair_data_from_context(ctx,row,payload_mode)
data = ctx.Value;
pair = build_pair_data(data.q,data.rbase_in_c,data.rbase_in_f, ...
    data.rimage_pairs,data.refine,data.pairs,data.opt,data.Lc,row, ...
    false,data.svd_opts,false,payload_mode);
end

function [Uf,Yf,Up,Yp,Cmap,Cmap_FU] = store_pair_payload( ...
    Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair,pairs,row,N_peanut,use_pair_map, ...
    store_full_pair_payload)
i = pairs(row,1);
p2 = pairs(row,2);

if store_full_pair_payload
    Uf{i,p2} = pair.Upf;
    Yf{i,p2} = pair.Ypf;
end

if N_peanut && ~use_pair_map
    if store_full_pair_payload
        Up{i,p2} = pair.DC;
        Yp{i,p2} = pair.YC;
    end
elseif N_peanut && use_pair_map
    Cmap{i,p2} = pair.Cmap;
end

if use_pair_map
    Cmap_FU{i,p2} = pair.Cmap_FU;
end
end

function group = init_pair_group()
group = struct('group_id',[],'sep',[],'q_pair',[],'rimage_canon',{{}}, ...
    'refine_canon',{{}},'Upf',[],'Ypf',[],'DC',[],'YC',[], ...
    'Cmap',[],'Cmap_proj',[],'Cmap_FU',[],'Lf_pair',[],'Lc_pair',[], ...
    'Upair_colloc',[],'Ucross_colloc',[],'Ecolloc',[],'rep_pair',[]);
end

function meta = build_pair_meta(q,pairs,nc,nf,opt)
total_pairs = size(pairs,1);
meta = repmat(struct('i',[],'j',[],'group_id',[],'sep',[], ...
    'mid',[],'rot',[],'phase_c',[],'phase_c_inv',[], ...
    'phase_f',[],'phase_f_inv',[], ...
    'Upair_colloc_actual',[],'Ucross_colloc_actual',[], ...
    'Ecolloc_actual',[]),total_pairs,1);

for row = 1:total_pairs
    i = pairs(row,1);
    j = pairs(row,2);
    delta = q(j)-q(i);
    sep = abs(delta);
    if sep == 0
        rot = 1;
    else
        rot = delta/sep;
    end

    meta(row).i = i;
    meta(row).j = j;
    meta(row).sep = sep;
    meta(row).mid = 0.5*(q(i)+q(j));
    meta(row).rot = rot;
    meta(row).phase_c = getUniformCircleRotationSpec(nc,rot,opt);
    meta(row).phase_c_inv = invertUniformCircleRotationSpec(meta(row).phase_c);
    meta(row).phase_f = getUniformCircleRotationSpec(nf,rot,opt);
    meta(row).phase_f_inv = invertUniformCircleRotationSpec(meta(row).phase_f);
end

end

function [group_id,group_sep,rep_rows] = group_pair_separations(meta,sep_tol)
sep = reshape([meta.sep],[],1);
[sep_sorted,order] = sort(sep);

group_id_sorted = zeros(size(sep_sorted));
group_sep = zeros(0,1);
rep_rows = zeros(0,1);

if isempty(sep_sorted)
    group_id = zeros(0,1);
    return
end

group_count = 1;
group_id_sorted(1) = group_count;
group_sep(1,1) = sep_sorted(1);
rep_rows(1,1) = order(1);

for k = 2:numel(sep_sorted)
    if abs(sep_sorted(k)-group_sep(group_count)) <= sep_tol
        group_id_sorted(k) = group_count;
    else
        group_count = group_count + 1;
        group_id_sorted(k) = group_count;
        group_sep(group_count,1) = sep_sorted(k);
        rep_rows(group_count,1) = order(k);
    end
end

group_id = zeros(size(sep));
group_id(order) = group_id_sorted;
end

function group = build_pair_group(group_id,row,group_sep,q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,store_full_pair_payload,debug,svd_opts,payload_mode)
pair = build_pair_data(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,row,debug,svd_opts,true,payload_mode);
group = init_pair_group();
group.group_id = group_id;
group.sep = group_sep;
group.q_pair = pair.q_pair;
group.rimage_canon = pair.rimage_canon;
group.refine_canon = pair.refine_canon;
group.Cmap = pair.Cmap;
group.Cmap_FU = pair.Cmap_FU;
group.rep_pair = pair.rep_pair;

if store_full_pair_payload
    group.Upf = pair.Upf;
    group.Ypf = pair.Ypf;
    group.DC = pair.DC;
    group.YC = pair.YC;
    group.Lf_pair = pair.Lf_pair;
    group.Lc_pair = pair.Lc_pair;
    group.Upair_colloc = pair.Upair_colloc;
    group.Ucross_colloc = pair.Ucross_colloc;
    group.Ecolloc = pair.Ecolloc;
end
end

function pair = build_pair_data(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,row,debug,svd_opts,use_canonical,payload_mode)
if nargin < 11 || isempty(svd_opts)
    svd_opts = struct();
end
if nargin < 12 || isempty(use_canonical)
    use_canonical = false;
end
if nargin < 13 || isempty(payload_mode)
    payload_mode = 'full';
end
store_full_pair_payload = strcmp(payload_mode,'full');

i = pairs(row,1);
j = pairs(row,2);
pair.rep_pair = [i j];

delta = q(j)-q(i);
sep = abs(delta);
if sep == 0
    rot = 1;
else
    rot = delta/sep;
end
mid = 0.5*(q(i)+q(j));

if use_canonical
    q_pair = [-sep/2; sep/2];
    rimage_i = map_points_to_canonical(rimage_pairs{i,j},mid,rot);
    rimage_j = map_points_to_canonical(rimage_pairs{j,i},mid,rot);
    refine_i = map_points_to_canonical(refine{i,j},mid,rot);
    refine_j = map_points_to_canonical(refine{j,i},mid,rot);
else
    q_pair = [q(i); q(j)];
    rimage_i = rimage_pairs{i,j};
    rimage_j = rimage_pairs{j,i};
    refine_i = refine{i,j};
    refine_j = refine{j,i};
end

pair.q_pair = q_pair;
pair.rimage_canon = {rimage_i; rimage_j};
pair.refine_canon = {refine_i; refine_j};
pair.rin_pair = [];
pair.Upf = [];
pair.Ypf = [];
pair.Lf_pair = [];
pair.Lc_pair = [];
pair.DC = [];
pair.YC = [];
pair.Cmap = [];
pair.Cmap_FU = [];
pair.Upair_colloc = [];
pair.Ucross_colloc = [];
pair.Ecolloc = [];


project_force = opt.project_force;

N_c = opt.N_c;
a_c = getOptField(opt,'a_c',1.2);
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';
rout_base_c = cos(tout_c)+1i*sin(tout_c);

nout = ceil(opt.a_f*opt.N_f);
t = linspace(0,2*pi,nout+1);
t = t(1:end-1)';
rout_base = cos(t)+1i*sin(t);

rin_1_f = q_pair(1)+rbase_in_f;
rin_2_f = q_pair(2)+rbase_in_f;
rout_f = [q_pair(1)+rout_base; refine_i; q_pair(2)+rout_base; refine_j];
rin_pair = [rin_1_f; rimage_i; rin_2_f; rimage_j];
if store_full_pair_payload
    pair.rin_pair = rin_pair;
end

svd_pair = svd_opts;
if logical(getOptField(svd_opts,'left_weight',false))
    row_weights_1 = getPeriodicCurveWeights([q_pair(1)+rout_base; refine_i],q_pair(1));
    row_weights_2 = getPeriodicCurveWeights([q_pair(2)+rout_base; refine_j],q_pair(2));
    svd_pair.row_weights = [row_weights_1; row_weights_2];
end

if debug
    figure(801);
    clf;
    plot(real(rin_1_f),imag(rin_1_f),'r.','MarkerSize',10);
    hold on;
    plot(real(rin_2_f),imag(rin_2_f),'b.','MarkerSize',10);
    plot(real(q_pair(1)+rout_base),imag(q_pair(1)+rout_base),'ro','MarkerSize',4);
    plot(real(q_pair(2)+rout_base),imag(q_pair(2)+rout_base),'bo','MarkerSize',4);
    plot(real(refine_i),imag(refine_i),'r+','MarkerSize',6);
    plot(real(refine_j),imag(refine_j),'b+','MarkerSize',6);
    if ~isempty(rimage_i)
        plot(real(rimage_i),imag(rimage_i),'ks','MarkerSize',5);
    end
    if ~isempty(rimage_j)
        plot(real(rimage_j),imag(rimage_j),'kd','MarkerSize',5);
    end
    axis equal;
    grid on;
    title(sprintf('getPairBasisStokes pair (%d,%d)',i,j), ...
        'Interpreter','none');
    drawnow;
end

if ~isempty(Lc)
    Lc_pair = getILpair(Lc);
else
    Lc_pair = [];
end

if project_force || opt.cmap
    Kf1 = getKmat2D(rin_pair(1:end/2),q_pair(1));
    Kf2 = getKmat2D(rin_pair(end/2+1:end),q_pair(2));
else
    Kf1 = [];
    Kf2 = [];
end

if project_force
    B1 = getKmat2D([q_pair(1)+rout_base; refine_i],q_pair(1));
    B2 = getKmat2D([q_pair(2)+rout_base; refine_j],q_pair(2));
    Lr_pair = getLrPair(B1,B2,Kf1,Kf2);
    Lf_pair = getLfPair(Kf1,Kf2);
else
    Lr_pair = [];
    Lf_pair = [];
end

[Uf_pair,Yf_pair] = getPairBlockStokes(rin_pair,rout_f,Lf_pair,Lr_pair,svd_pair);
Npair = evaluateCoarseOnPair(q_pair,rbase_in_c,rout_f);
Upf = -Uf_pair'*Npair; %

DC = [];
YC = [];
if opt.N_peanut
    rout_peanut = createPeanut(q_pair(1),q_pair(2),opt.N_peanut,0);
    rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
    [DC,YC] = getPeanutBlockStokes(rin_pair_c,rin_pair,rout_peanut,Lc_pair,Lf_pair,svd_opts);
    if opt.cmap
        pair.Cmap = -YC*(DC*Yf_pair*(Uf_pair'*Npair));
    end
else
    rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
end


if opt.cmap
    Kft_pair = getKftPair(Kf1,Kf2);
    pair.Cmap_FU = -Kft_pair*Yf_pair*(Uf_pair'*Npair);
end

if store_full_pair_payload
    pair.Upf = Upf;
    pair.Ypf = Yf_pair;
    pair.Lf_pair = Lf_pair;
    pair.Lc_pair = Lc_pair;
    pair.DC = DC;
    pair.YC = YC;
end

if debug
    run_pair_lsq_debug_test(i,j,row,q_pair,rbase_in_c,rin_pair,Upf,Yf_pair, ...
        Lf_pair,Kf1,Kf2,Lc_pair,DC,YC,rout_f,Npair,opt);
end

if store_full_pair_payload && project_force
    rout_pair_c = [q_pair(1)+rout_base_c; q_pair(2)+rout_base_c];
    Ppair = eye(size(Lf_pair)) - Lf_pair;
    pair.Upair_colloc = stokSLPmat(rin_pair,rout_pair_c,1) * Ppair * Yf_pair * Upf;
    pair.Ecolloc = stokSLPmat(rin_pair_c,rout_pair_c,1);
    pair.Ucross_colloc = build_cross_pair_velocity_map(pair.Ecolloc,N_c,numel(rout_base_c));
end
end

function run_pair_lsq_debug_test(i,j,row,q_pair,rbase_in_c,rin_pair,Upf,Yf_pair, ...
    Lf_pair,Kf1,Kf2,Lc_pair,DC,YC,rout_pair,Npair_colloc,opt)
% Debug-only check of the two LS problems used in getPairBasisStokes.
%
% LS1 right-hand side is denoted r_pair in this file:
%   r_pair = -Npair * [lambda_1x; lambda_2x; lambda_1y; lambda_2y]
% where lambda_1 and lambda_2 are coarse source vectors on each particle.
%
% LS2 uses Ntot:
%   Ntot = Nf - Nf*Lf_pair  (mobility)
%   Ntot = Nf               (resistance)

N_c = opt.N_c;
[lambda_1,lambda_2] = build_pair_debug_rhs(N_c,row,opt);
rhs_pair = [lambda_1(1:N_c); lambda_2(1:N_c); ...
    lambda_1(N_c+1:end); lambda_2(N_c+1:end)];

pair_mapped = Upf*rhs_pair;
beta = Yf_pair*pair_mapped;

nout_test = 2*length(rout_pair);
ntest_half = nout_test/2;
rout_pair_test = build_gap_refined_pair_boundary_nodes(q_pair,ntest_half);

S_test = stokSLPmat(rin_pair,rout_pair_test,1);
if isempty(Lf_pair)
    A_test = S_test;
else
    B1_test = getKmat2D(rout_pair_test(1:ntest_half),q_pair(1));
    B2_test = getKmat2D(rout_pair_test(ntest_half+1:end),q_pair(2));
    Lr_test = getLrPair(B1_test,B2_test,Kf1,Kf2);
    A_test = S_test-S_test*Lf_pair+Lr_test;
end

Npair_test = evaluateCoarseOnPair(q_pair,rbase_in_c,rout_pair_test);
r_pair_test = -Npair_test*rhs_pair;
ls1_res = A_test*beta-r_pair_test;
[ls1_abs_inf,ls1_abs_pw,ls1_rel,ls1_rel_pw] = pairwise_error_metrics(ls1_res,r_pair_test);

S_colloc = stokSLPmat(rin_pair,rout_pair,1);
if isempty(Lf_pair)
    A_colloc = S_colloc;
else
    B1_colloc = getKmat2D(rout_pair(1:end/2),q_pair(1));
    B2_colloc = getKmat2D(rout_pair(end/2+1:end),q_pair(2));
    Lr_colloc = getLrPair(B1_colloc,B2_colloc,Kf1,Kf2);
    A_colloc = S_colloc-S_colloc*Lf_pair+Lr_colloc;
end
r_pair_colloc = -Npair_colloc*rhs_pair;
ls1_res_colloc = A_colloc*beta-r_pair_colloc;
[ls1_abs_inf_colloc,ls1_abs_pw_colloc,ls1_rel_colloc,ls1_rel_pw_colloc] = ...
    pairwise_error_metrics(ls1_res_colloc,r_pair_colloc);

fig1 = 8400 + 2*(row-1) + 1;
plot_abs_rel_error_compare(fig1,ls1_abs_pw,ls1_abs_pw_colloc, ...
    ls1_rel_pw,ls1_rel_pw_colloc,'Boundary node number', ...
    sprintf(['Pair (%d,%d) LSQ for fine sources: abs_{inf} off/coll = %.3e / %.3e, ', ...
    'rel_{inf} off/coll = %.3e / %.3e'], ...
    i,j,ls1_abs_inf,ls1_abs_inf_colloc,ls1_rel,ls1_rel_colloc));

ls2_rel = NaN;
ls2_rel_colloc = NaN;
ls2_abs_inf = NaN;
ls2_abs_inf_colloc = NaN;
if ~isempty(DC) && ~isempty(YC) && (getOptField(opt,'N_peanut',0) > 0)
    rout_peanut_colloc = createPeanut(q_pair(1),q_pair(2),opt.N_peanut,0);
    n_peanut_test = 2*length(rout_peanut_colloc);
    rout_peanut_test = createPeanut(q_pair(1),q_pair(2),n_peanut_test,0);

    rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
    Nf_test = stokSLPmat(rin_pair,rout_peanut_test,1);
    if isempty(Lf_pair)
        Ntot_test = Nf_test;
    else
        Ntot_test = Nf_test-Nf_test*Lf_pair;
    end

    Npeanut_test = stokSLPmat(rin_pair_c,rout_peanut_test,1);
    if ~isempty(Lc_pair)
        Npeanut_test = Npeanut_test*Lc_pair;
    end

    tau_peanut = YC*(DC*beta);
    rhs2 = Ntot_test*beta;
    ls2_res = Npeanut_test*tau_peanut-rhs2;
    [ls2_abs_inf,ls2_abs_pw,ls2_rel,ls2_rel_pw] = pairwise_error_metrics(ls2_res,rhs2);

    Nf_colloc = stokSLPmat(rin_pair,rout_peanut_colloc,1);
    if isempty(Lf_pair)
        Ntot_colloc = Nf_colloc;
    else
        Ntot_colloc = Nf_colloc-Nf_colloc*Lf_pair;
    end
    Npeanut_colloc = stokSLPmat(rin_pair_c,rout_peanut_colloc,1);
    if ~isempty(Lc_pair)
        Npeanut_colloc = Npeanut_colloc*Lc_pair;
    end
    rhs2_colloc = Ntot_colloc*beta;
    ls2_res_colloc = Npeanut_colloc*tau_peanut-rhs2_colloc;
    [ls2_abs_inf_colloc,ls2_abs_pw_colloc,ls2_rel_colloc,ls2_rel_pw_colloc] = ...
        pairwise_error_metrics(ls2_res_colloc,rhs2_colloc);

    fig2 = fig1 + 1;
    plot_abs_rel_error_compare(fig2,ls2_abs_pw,ls2_abs_pw_colloc, ...
        ls2_rel_pw,ls2_rel_pw_colloc,'boundary node number', ...
        sprintf(['Pair (%d,%d) peanut LSQ problem: abs_{inf} off/coll = %.3e / %.3e, ', ...
        'rel_{inf} off/coll = %.3e / %.3e'], ...
        i,j,ls2_abs_inf,ls2_abs_inf_colloc,ls2_rel,ls2_rel_colloc));
end

fprintf(['getPairBasisStokes debug pair (%d,%d): ', ...
    'LSQ fine absinf off/coll = %.3e / %.3e, LS1 relinf off/coll = %.3e / %.3e, ', ...
    'LSQ peanut absinf off/coll = %.3e / %.3e, LS2 relinf off/coll = %.3e / %.3e\n'], ...
    i,j,ls1_abs_inf,ls1_abs_inf_colloc,ls1_rel,ls1_rel_colloc, ...
    ls2_abs_inf,ls2_abs_inf_colloc,ls2_rel,ls2_rel_colloc);
end

function [lambda_1,lambda_2] = build_pair_debug_rhs(N_c,row,opt)
has_lam1 = isfield(opt,'pair_basis_lambda_1') && ~isempty(opt.pair_basis_lambda_1);
has_lam2 = isfield(opt,'pair_basis_lambda_2') && ~isempty(opt.pair_basis_lambda_2);
if has_lam1 || has_lam2
    if ~(has_lam1 && has_lam2)
        error(['pair_basis_lambda_1 and pair_basis_lambda_2 must be set ', ...
            'together when debug is enabled.']);
    end
    lambda_1 = opt.pair_basis_lambda_1(:);
    lambda_2 = opt.pair_basis_lambda_2(:);
    if numel(lambda_1) ~= 2*N_c || numel(lambda_2) ~= 2*N_c
        error('pair_basis_lambda_1 and pair_basis_lambda_2 must each have length 2*N_c.');
    end
    return
end

t = linspace(0,2*pi,N_c+1)';
t = t(1:end-1);
phi = 0.37*row;

lambda1_x = cos(t+phi) + 0.35*sin(2*t-0.2) + 0.15*cos(3*t+0.7);
lambda1_y = 0.8*lambda1_x - 0.3*sin(3*t+phi);
lambda2_x = sin(t-0.3*phi) - 0.25*cos(2*t+0.4) + 0.12*sin(4*t-0.9);
lambda2_y = -0.6*lambda2_x + 0.2*cos(2*t-phi);

lambda_1 = [lambda1_x; lambda1_y];
lambda_2 = [lambda2_x; lambda2_y];
end

function rout_pair = build_gap_refined_pair_boundary_nodes(q_pair,nout)
% Build off-collocation check nodes with:
% 1) a uniform set over the whole boundary, and
% 2) an extra local cluster near the pair gap using a Mobius angle map.
% Gap point on each particle is the point facing the other center.
delta = q_pair(2)-q_pair(1);
if abs(delta) == 0
    gap_theta_1 = 0;
else
    gap_theta_1 = angle(delta);
end
gap_theta_2 = mod(gap_theta_1+pi,2*pi);

% Split between globally uniform nodes and locally clustered nodes.
n_uniform = max(8,round(0.45*nout));
n_uniform = min(n_uniform,nout);
n_cluster = nout - n_uniform;

tu = linspace(0,2*pi,n_uniform+1)';
tu = tu(1:end-1) + pi/n_uniform;
tu = mod(tu,2*pi);

tc = zeros(0,1);
if n_cluster > 0
    % Build a Mobius-mapped angle pool and keep the nodes closest to 0,
    % i.e. closest to the local gap direction.
    r_mob = 0.88;
    npool = max(8*n_cluster,128);
    tseed = linspace(0,2*pi,npool+1)';
    tseed = tseed(1:end-1);
    tmob = sort(mobius_angle_map(tseed,r_mob));
    dgap = min(abs(tmob),2*pi-abs(tmob));
    [~,ord] = sort(dgap,'ascend');
    tc = tmob(ord(1:n_cluster));
end

t_local = sort(mod([tu; tc],2*pi));

t1 = mod(gap_theta_1 + t_local,2*pi);
t2 = mod(gap_theta_2 + t_local,2*pi);
rout_1 = q_pair(1) + exp(1i*t1);
rout_2 = q_pair(2) + exp(1i*t2);
rout_pair = [rout_1; rout_2];
end

function theta = mobius_angle_map(t,r)
% Same Mobius map form used in geometry/pair_clusters_ellipse.m.
z = exp(1i*t);
w = (z-r)./(1-r*z);
theta = mod(angle(w),2*pi);
end

function [abs_inf,abs_pw,rel_inf,rel_pw] = pairwise_error_metrics(err_vec,ref_vec)
n = numel(err_vec)/2;
err_x = err_vec(1:n);
err_y = err_vec(n+1:end);
ref_x = ref_vec(1:n);
ref_y = ref_vec(n+1:end);

abs_pw = hypot(err_x,err_y);
abs_inf = max(abs_pw);
ref_mag = hypot(ref_x,ref_y);
ref_scale = max(ref_mag);
den = max(ref_mag,1e-14*max(1,ref_scale));
rel_pw = abs_pw./den;
rel_inf = norm(err_vec,inf)/max(1,norm(ref_vec,inf));
end

function plot_abs_rel_error_compare(fig_id,abs_off,abs_coll,rel_off,rel_coll,xlab,tstr)
figure(fig_id);
clf;
subplot(2,1,1);
semilogy(abs_off + eps,'b-','LineWidth',1.2);
hold on;
semilogy(abs_coll + eps,'k--','LineWidth',1.2);
grid on;
ylabel('Absolute error');
legend('new nodes','collocation','Location','best');
title(tstr,'Interpreter','none');
axis tight;

subplot(2,1,2);
semilogy(rel_off + eps,'r-','LineWidth',1.2);
hold on;
semilogy(rel_coll + eps,'m--','LineWidth',1.2);
grid on;
xlabel(xlab);
ylabel('Relative error');
legend('new nodes','collocation','Location','best');
axis tight;

end

function pair_cache = populate_actual_dense_pair_fields(pair_cache,q,rbase_in_c,rbase_in_f, ...
    rout_base_c,rimage_pairs,pairs)
total_pairs = size(pairs,1);
for row = 1:total_pairs
    [Upair_actual,Ucross_actual,E_actual] = build_actual_pair_dense( ...
        q,rbase_in_c,rbase_in_f,rout_base_c,rimage_pairs,pairs,row);
    pair_cache.meta(row).Upair_colloc_actual = Upair_actual;
    pair_cache.meta(row).Ucross_colloc_actual = Ucross_actual;
    pair_cache.meta(row).Ecolloc_actual = E_actual;
end
end

function [Upair_actual,Ucross_actual,E_actual] = build_actual_pair_dense( ...
    q,rbase_in_c,rbase_in_f,rout_base_c,rimage_pairs,pairs,row)
i = pairs(row,1);
j = pairs(row,2);

rin_pair_f = [q(i)+rbase_in_f; rimage_pairs{i,j}; q(j)+rbase_in_f; rimage_pairs{j,i}];
rin_pair_c = [q(i)+rbase_in_c; q(j)+rbase_in_c];
rout_pair = [q(i)+rout_base_c; q(j)+rout_base_c];

Upair_actual = stokSLPmat(rin_pair_f,rout_pair,1);
E_actual = stokSLPmat(rin_pair_c,rout_pair,1);
Ucross_actual = build_cross_pair_velocity_map(E_actual,numel(rbase_in_c),numel(rout_base_c));
end

function Ucross = build_cross_pair_velocity_map(Epair,N_src,N_tgt)
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

function z_canon = map_points_to_canonical(z,mid,rot)
z = z(:);
if isempty(z)
    z_canon = zeros(0,1);
    return
end
z_canon = conj(rot)*(z-mid);
end
