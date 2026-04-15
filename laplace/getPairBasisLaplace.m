function [Uf,Yf,Up,Yp,Cmap,Cmap_QV,pair_cache] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f,rout_base_c,rimage_vec,refine,pairs,opt)
%GETPAIRBASISLAPLACE Build pair-basis pseudoinverse factors for Laplace.
%
% Syntax:
%   [Uf,Yf,Up,Yp,Cmap,Cmap_QV,pair_cache] = ...
%       getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f,rout_base_c, ...
%       rimage_vec,refine,pairs,opt)
%
% See also: getPairBlockLaplace, getPeanutBlockLaplace, ...
%   evaluateCoarseOnPairLaplace, getPairTransformationLaplace.
%
% Anna Broms, Mar 2026

q = q(:);
P = numel(q);
N_f = opt.N_f;
N_peanut = opt.N_peanut;
R = opt.rad;
if numel(R) > 1
    R = R(1);
end
reuse_pair_basis = opt.reuse_pair_basis_by_sep;
parallel_precomp = opt.parallel_precomp;
check_rotations = opt.check_rotations;
shared_sep_tol = opt.shared_sep_tol;
project_charge = opt.project_charge;
show_counter = opt.show_counter;
use_pair_map = opt.cmap;

if reuse_pair_basis
    Uf = [];
    Yf = [];
    Up = [];
    Yp = [];
    Cmap = [];
    Cmap_QV = [];
else
    Uf = cell(P);
    Yf = cell(P);

    if N_peanut
        if use_pair_map
            Up = [];
            Yp = [];
            Cmap = cell(P);
            Cmap_QV = cell(P);
        else
            Up = cell(P);
            Yp = cell(P);
            Cmap = [];
            Cmap_QV = [];
        end
    else
        Up = [];
        Yp = [];
        Cmap = [];
        Cmap_QV = [];
    end
end

pair_cache = init_pair_cache();
pair_cache.enabled = reuse_pair_basis;
pair_cache.check_rotations = check_rotations;
pair_cache.shared_sep_tol = shared_sep_tol;
pair_cache.rout_base_f = rout_base_f(:);
pair_cache.stats.requested_parallel = parallel_precomp;

if isempty(pairs)
    return
end

pair_cache.meta = build_pair_meta(q,pairs,numel(rbase_in_c),numel(rbase_in_f),opt);
pair_cache.stats.n_pairs = size(pairs,1);

if ~reuse_pair_basis
    total_pairs = size(pairs,1);
    pair_cache.stats.branch = 'per_pair';
    use_parallel_pairs = parallel_precomp && total_pairs > 1;
    [use_parallel_pairs,pool_size] = resolve_parallel_precomp(use_parallel_pairs, ...
        'getPairBasisLaplace');
    pair_cache.stats.used_parallel = use_parallel_pairs;
    pair_cache.stats.pool_size = pool_size;

    if use_parallel_pairs
        if show_counter
            fprintf('getPairBasisLaplace: parallel pair build for %d pairs\n', total_pairs);
        end
        pair_entries = cell(total_pairs,1);
        parfor ii = 1:total_pairs
            pair_entries{ii} = build_pair_entry(ii,q,rbase_in_c,rbase_in_f,rimage_vec, ...
                refine,pairs,opt,pair_cache.rout_base_f,project_charge, ...
                pair_cache.meta(ii).sep);
        end
        pair_cache.stats.pool_size = get_parallel_pool_size();
        for ii = 1:total_pairs
            [Uf,Yf,Up,Yp,Cmap,Cmap_QV] = assign_pair_entry_outputs( ...
                Uf,Yf,Up,Yp,Cmap,Cmap_QV,pairs,ii,pair_entries{ii},N_peanut,use_pair_map);
        end
        if show_counter
            fprintf('getPairBasisLaplace: finished parallel pair build for %d pairs\n', ...
                total_pairs);
        end
    else
        processed_pairs = 0;
        for ii = 1:total_pairs
            pair_entry = build_pair_entry(ii,q,rbase_in_c,rbase_in_f,rimage_vec, ...
                refine,pairs,opt,pair_cache.rout_base_f,project_charge, ...
                pair_cache.meta(ii).sep);
            [Uf,Yf,Up,Yp,Cmap,Cmap_QV] = assign_pair_entry_outputs( ...
                Uf,Yf,Up,Yp,Cmap,Cmap_QV,pairs,ii,pair_entry,N_peanut,use_pair_map);

            processed_pairs = processed_pairs + 1;
            if show_counter
                i = pairs(ii,1);
                p2 = pairs(ii,2);
                fprintf('getPairBasisLaplace: processed pair %d/%d (%d,%d)\n', ...
                    processed_pairs,total_pairs,i,p2);
            end
        end
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
    'getPairBasisLaplace');
pair_cache.stats.used_parallel = use_parallel_groups;
pair_cache.stats.pool_size = pool_size;

if use_parallel_groups
    if show_counter
        fprintf(['getPairBasisLaplace: parallel canonical build for %d groups ', ...
            'covering %d pairs\n'], n_groups,total_pairs);
    end
    groups = repmat(init_pair_group(),n_groups,1);
    parfor gg = 1:n_groups
        groups(gg) = build_pair_group(gg,rep_rows(gg),group_sep(gg), ...
            q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,pair_cache.rout_base_f, ...
            project_charge);
    end
    pair_cache.groups = groups;
    pair_cache.stats.pool_size = get_parallel_pool_size();
else
    pair_cache.groups = repmat(init_pair_group(),n_groups,1);
    covered_pairs = 0;
    for gg = 1:n_groups
        pair_cache.groups(gg) = build_pair_group(gg,rep_rows(gg),group_sep(gg), ...
            q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,pair_cache.rout_base_f, ...
            project_charge);

        if show_counter
            rep_i = pairs(rep_rows(gg),1);
            rep_j = pairs(rep_rows(gg),2);
            n_rep = sum(group_id == gg);
            covered_pairs = covered_pairs + n_rep;
            fprintf(['getPairBasisLaplace: processed canonical group %d/%d ', ...
                'from pair (%d,%d), pp_sep = %.3g, covers %d pairs -> %d/%d pairs covered\n'], ...
                gg,n_groups,rep_i,rep_j,group_sep(gg)-2*R,n_rep,covered_pairs,total_pairs);
        end
    end
end

if check_rotations
    % Keep a per-pair copy in the solve geometry so the debug path can
    % compare against the canonical cached group later.
    pair_cache.check_pairs = repmat(init_pair_group(),size(pairs,1),1);
    for row = 1:size(pairs,1)
        gid = group_id(row);
        pair_cache.check_pairs(row) = build_pair_group(gid,row,pair_cache.meta(row).sep, ...
            q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,pair_cache.rout_base_f, ...
            project_charge,false);
    end
end

for row = 1:size(pairs,1)
    gid = group_id(row);
    pair_cache.meta(row).group_id = gid;
    pair_cache.meta(row).sep = group_sep(gid);
    if use_pair_map
        [Ucross_actual,Ec_actual,Lr_actual] = build_actual_pair_colloc_factors( ...
            pair_cache.meta(row),q,rbase_in_c,rout_base_c);
        pair_cache.meta(row).Ucross_colloc_actual = Ucross_actual;
        pair_cache.meta(row).Ec_colloc_actual = Ec_actual;
        pair_cache.meta(row).Lr_colloc_actual = Lr_actual;
    end
end

if show_counter
    fprintf('getPairBasisLaplace: reused %d canonical groups for %d pairs\n', ...
        n_groups,size(pairs,1));
end

end

function pair_cache = init_pair_cache()
pair_cache = struct();
pair_cache.enabled = false;
pair_cache.check_rotations = false;
pair_cache.shared_sep_tol = [];
pair_cache.rout_base_f = [];
pair_cache.meta = repmat(struct('i',[],'j',[],'group_id',[],'sep',[], ...
    'mid',[],'rot',[],'phase_c',[],'phase_c_inv',[], ...
    'phase_f',[],'phase_f_inv',[], ...
    'Ucross_colloc_actual',[],'Ec_colloc_actual',[], ...
    'Lr_colloc_actual',[]),0,1);
pair_cache.groups = repmat(init_pair_group(),0,1);
pair_cache.check_pairs = repmat(init_pair_group(),0,1);
pair_cache.group_id = zeros(0,1);
pair_cache.group_sep = zeros(0,1);
pair_cache.representative_rows = zeros(0,1);
pair_cache.n_groups = 0;
pair_cache.stats = init_pair_stats();
end

function stats = init_pair_stats()
stats = struct('requested_parallel',false,'used_parallel',false, ...
    'branch','','n_pairs',0,'n_groups',0,'pool_size',0);
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

pool_size = get_parallel_pool_size();
end

function pool_size = get_parallel_pool_size()
pool_size = 0;
pool = gcp('nocreate');
if ~isempty(pool)
    pool_size = pool.NumWorkers;
end
end

function group = init_pair_group()
group = struct('group_id',[],'sep',[],'q_pair',[],'rimage_canon',{{}}, ...
    'refine_canon',{{}},'Upf',[],'Ypf',[],'DC',[],'YC',[],'Cmap',[], ...
    'Cmap_QV',[], ...
    'nout_f',[],'nsrc_f',[],'ntar_f',[],'rep_pair',[]);
end

function meta = build_pair_meta(q,pairs,nc,nf,opt)
total_pairs = size(pairs,1);
meta = repmat(struct('i',[],'j',[],'group_id',[],'sep',[], ...
    'mid',[],'rot',[],'phase_c',[],'phase_c_inv',[], ...
    'phase_f',[],'phase_f_inv',[], ...
    'Ucross_colloc_actual',[],'Ec_colloc_actual',[], ...
    'Lr_colloc_actual',[]),total_pairs,1);

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

function pair_entry = build_pair_entry(row,q,rbase_in_c,rbase_in_f,rimage_vec, ...
    refine,pairs,opt,rout_base_f,project_charge,sep)
pair_entry = build_pair_group([],row,sep,q,rbase_in_c,rbase_in_f,rimage_vec, ...
    refine,pairs,opt,rout_base_f,project_charge,false);
end

function [Uf,Yf,Up,Yp,Cmap,Cmap_QV] = assign_pair_entry_outputs( ...
    Uf,Yf,Up,Yp,Cmap,Cmap_QV,pairs,row,pair_entry,N_peanut,use_pair_map)
i = pairs(row,1);
j = pairs(row,2);
Uf{i,j} = pair_entry.Upf;
Yf{i,j} = pair_entry.Ypf;

if N_peanut
    if use_pair_map
        Cmap{i,j} = pair_entry.Cmap;
        Cmap_QV{i,j} = pair_entry.Cmap_QV;
    else
        Up{i,j} = pair_entry.DC;
        Yp{i,j} = pair_entry.YC;
    end
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

function group = build_pair_group(group_id,row,sep,q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,rout_base_f,project_charge,use_canonical)
if nargin < 13
    use_canonical = true;
end

i = pairs(row,1);
j = pairs(row,2);
mid = 0.5*(q(i)+q(j));
rot = (q(j)-q(i))/abs(q(j)-q(i));
R = opt.rad;

if use_canonical
    q_pair = [-sep/2; sep/2];
    rimage_i = map_points_to_canonical(rimage_vec{i,j},mid,rot);
    rimage_j = map_points_to_canonical(rimage_vec{j,i},mid,rot);
    refine_i = map_points_to_canonical(refine{i,j},mid,rot);
    refine_j = map_points_to_canonical(refine{j,i},mid,rot);
else
    q_pair = [q(i); q(j)];
    rimage_i = rimage_vec{i,j};
    rimage_j = rimage_vec{j,i};
    refine_i = refine{i,j};
    refine_j = refine{j,i};
end

rout_f = [q_pair(1)+rout_base_f; refine_i; q_pair(2)+rout_base_f; refine_j];
rin_pair_f = [q_pair(1)+rbase_in_f; rimage_i; q_pair(2)+rbase_in_f; rimage_j];
nsrc_f_i = numel(rbase_in_f) + numel(rimage_i);
nsrc_f_j = numel(rbase_in_f) + numel(rimage_j);
ntar_i = numel(rout_base_f) + numel(refine_i);
ntar_j = numel(rout_base_f) + numel(refine_j);

proj_pair = struct('project_charge',project_charge, ...
    'nsrc',[nsrc_f_i nsrc_f_j], ...
    'ntar',[ntar_i ntar_j]);
[Uf_pair,Yf_pair] = getPairBlockLaplace(rin_pair_f,rout_f,proj_pair);

Npair = evaluateCoarseOnPairLaplace(q_pair,rbase_in_c,rout_f);
Upf = -Uf_pair'*Npair;
Ypf = Yf_pair;

DC = [];
YC = [];
C = [];
QV = [];

if opt.N_peanut
    rout_peanut = createPeanut(q_pair(1),q_pair(2),opt.N_peanut,0,R);
    rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
    proj_peanut = struct('project_charge',project_charge, ...
        'nsrc_c',[numel(rbase_in_c) numel(rbase_in_c)], ...
        'nsrc_f',[nsrc_f_i nsrc_f_j]);
    [DC,YC] = getPeanutBlockLaplace(rin_pair_c,rin_pair_f,rout_peanut,proj_peanut);

    if opt.cmap
        C = -YC*(DC*Yf_pair*(Uf_pair'*Npair));
        QV = getPairChargeSumMap(nsrc_f_i,nsrc_f_j)*Ypf*Upf;
        if opt.compress_cmap
            [U,S,V] = svd(C);
            S = diag(S);
            ra = sum(S>max(S)*opt.cmap_tol);
            fprintf('Rank of coarse-coarse map is chosen to %u\n',ra);
            S = S(1:ra);
            C = U(:,1:ra)*diag(S)*V(:,1:ra)';
        end
    end
end

group = init_pair_group();
group.group_id = group_id;
group.sep = sep;
group.q_pair = q_pair;
group.rimage_canon = {rimage_i; rimage_j};
group.refine_canon = {refine_i; refine_j};
group.Upf = Upf;
group.Ypf = Ypf;
group.DC = DC;
group.YC = YC;
group.Cmap = C;
group.Cmap_QV = QV;
group.nout_f = numel(rout_base_f);
group.nsrc_f = [nsrc_f_i nsrc_f_j];
group.ntar_f = [ntar_i ntar_j];
group.rep_pair = [i j];
end

function Kq_pair = getPairChargeSumMap(n1,n2)
Kq_pair = [ones(1,n1) zeros(1,n2); zeros(1,n1) ones(1,n2)];
end

function z_canon = map_points_to_canonical(z,mid,rot)
z = z(:);
if isempty(z)
    z_canon = zeros(0,1);
    return
end

z_canon = conj(rot)*(z-mid);
end

function [Ucross,Ec,Lr] = build_actual_pair_colloc_factors(meta,q,rbase_in_c,rout_base_c)

i = meta.i;
j = meta.j;
ntar_body = numel(rout_base_c);

rout_pair = [q(i)+rout_base_c; q(j)+rout_base_c];
rin_pair_c = [q(i)+rbase_in_c; q(j)+rbase_in_c];

Ec = lapSLPmat(rin_pair_c,rout_pair);
Ucross = -evaluateCoarseOnPairLaplace([q(i); q(j)],rbase_in_c,rout_pair);
Lr = build_pair_target_closure(ntar_body,ntar_body);
end

function Lr = build_pair_target_closure(m1,m2)

Lr = [ones(m1,1) zeros(m1,1); zeros(m2,1) ones(m2,1)];
end
