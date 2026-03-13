function [Uf,Yf,Up,Yp,Cmap,pair_cache] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f,rimage_vec,refine,pairs,opt)
%GETPAIRBASISLAPLACE Build pair-basis pseudoinverse factors for Laplace.
%
% Syntax:
%   [Uf,Yf,Up,Yp,Cmap,pair_cache] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f,rimage_vec,refine,pairs,opt)
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
reuse_pair_basis = logical(getOptField(opt,'reuse_pair_basis_by_sep',false));
shared_sep_tol = getOptField(opt,'shared_sep_tol',1e-12*max(1,R));
precomp = getOptField(opt,'precomp',true);

if isfield(opt,'project_charge') && ~isempty(opt.project_charge)
    project_charge = logical(opt.project_charge);
else
    project_charge = false;
end

if isfield(opt,'show_counter') && ~isempty(opt.show_counter)
    show_counter = logical(opt.show_counter);
else
    show_counter = false;
end

if reuse_pair_basis
    Uf = [];
    Yf = [];
    Up = [];
    Yp = [];
    Cmap = [];
else
    Uf = cell(P);
    Yf = cell(P);

    if N_peanut
        Up = cell(P);
        Yp = cell(P);
        Cmap = cell(P);
    else
        Up = [];
        Yp = [];
        Cmap = [];
    end
end

pair_cache = init_pair_cache();
pair_cache.enabled = reuse_pair_basis;
pair_cache.shared_sep_tol = shared_sep_tol;
pair_cache.rout_base_f = rout_base_f(:);

if isempty(pairs)
    return
end

pair_cache.meta = build_pair_meta(q,pairs);

if ~reuse_pair_basis
    total_pairs = size(pairs,1);
    processed_pairs = 0;

    for ii = 1:total_pairs
        i = pairs(ii,1);
        p2 = pairs(ii,2);

        fine_1 = refine{i,p2};
        fine_2 = refine{p2,i};
        rout_f = [q(i)+pair_cache.rout_base_f; fine_1; q(p2)+pair_cache.rout_base_f; fine_2];

        rimage_i = rimage_vec{i,p2};
        rimage_p2 = rimage_vec{p2,i};
        rin_pair_f = [q(i)+rbase_in_f; rimage_i; q(p2)+rbase_in_f; rimage_p2];
        nsrc_f_i = N_f + numel(rimage_i);
        nsrc_f_p2 = N_f + numel(rimage_p2);
        ntar_i = numel(pair_cache.rout_base_f) + numel(fine_1);
        ntar_p2 = numel(pair_cache.rout_base_f) + numel(fine_2);
        proj_pair = struct('project_charge',project_charge, ...
            'nsrc',[nsrc_f_i nsrc_f_p2], ...
            'ntar',[ntar_i ntar_p2]);

        [Uf_pair,Yf_pair] = getPairBlockLaplace(rin_pair_f,rout_f,proj_pair);

        Npair = evaluateCoarseOnPairLaplace([q(i);q(p2)],rbase_in_c,rout_f);
        if precomp
            Uf{i,p2} = -Uf_pair'*Npair;
        else
            Uf{i,p2} = Uf_pair';
        end
        Yf{i,p2} = Yf_pair;

        if N_peanut
            rout_peanut = createPeanut(q(i),q(p2),N_peanut,0,R);
            rin_pair_c = [q(i)+rbase_in_c; q(p2)+rbase_in_c];
            proj_peanut = struct('project_charge',project_charge, ...
                'nsrc_c',[numel(rbase_in_c) numel(rbase_in_c)], ...
                'nsrc_f',[nsrc_f_i nsrc_f_p2]);
            [DC,YC] = getPeanutBlockLaplace(rin_pair_c,rin_pair_f,rout_peanut,proj_peanut);

            if isfield(opt,'cmap') && opt.cmap
                C = -YC*(DC*Yf_pair*(Uf_pair'*Npair));
                if opt.compress_cmap
                    [U,S,V] = svd(C);
                    S = diag(S);
                    ra = sum(S>max(S)*opt.cmap_tol);
                    fprintf('Rank of coarse-coarse map is chosen to %u\n',ra);
                    S = S(1:ra);
                    C = U(:,1:ra)*diag(S)*V(:,1:ra)';
                end
                Cmap{i,p2} = C;
            else
                Up{i,p2} = DC;
                Yp{i,p2} = YC;
            end
        end

        processed_pairs = processed_pairs + 1;
        if show_counter
            fprintf('getPairBasisLaplace: processed pair %d/%d (%d,%d)\n', ...
                processed_pairs,total_pairs,i,p2);
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
pair_cache.groups = repmat(init_pair_group(),n_groups,1);

for gg = 1:n_groups
    pair_cache.groups(gg) = build_pair_group(gg,rep_rows(gg),group_sep(gg), ...
        q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,pair_cache.rout_base_f,project_charge,precomp);
end

for row = 1:size(pairs,1)
    gid = group_id(row);
    pair_cache.meta(row).group_id = gid;
    pair_cache.meta(row).sep = group_sep(gid);
end

if show_counter
    fprintf('getPairBasisLaplace: reused %d canonical groups for %d pairs\n', ...
        n_groups,size(pairs,1));
end

end

function pair_cache = init_pair_cache()
pair_cache = struct();
pair_cache.enabled = false;
pair_cache.shared_sep_tol = [];
pair_cache.rout_base_f = [];
pair_cache.meta = repmat(struct('i',[],'j',[],'group_id',[],'sep',[], ...
    'mid',[],'rot',[]),0,1);
pair_cache.groups = repmat(init_pair_group(),0,1);
pair_cache.group_id = zeros(0,1);
pair_cache.group_sep = zeros(0,1);
pair_cache.representative_rows = zeros(0,1);
pair_cache.n_groups = 0;
end

function group = init_pair_group()
group = struct('group_id',[],'sep',[],'q_pair',[],'rimage_canon',{{}}, ...
    'refine_canon',{{}},'Upf',[],'Ypf',[],'DC',[],'YC',[],'Cmap',[], ...
    'nout_f',[],'nsrc_f',[],'ntar_f',[],'rep_pair',[]);
end

function meta = build_pair_meta(q,pairs)
total_pairs = size(pairs,1);
meta = repmat(struct('i',[],'j',[],'group_id',[],'sep',[], ...
    'mid',[],'rot',[]),total_pairs,1);

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

function group = build_pair_group(group_id,row,sep,q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,rout_base_f,project_charge,precomp)
i = pairs(row,1);
j = pairs(row,2);
mid = 0.5*(q(i)+q(j));
rot = (q(j)-q(i))/abs(q(j)-q(i));
R = opt.rad;

q_pair = [-sep/2; sep/2];
rimage_i = map_points_to_canonical(rimage_vec{i,j},mid,rot);
rimage_j = map_points_to_canonical(rimage_vec{j,i},mid,rot);
refine_i = map_points_to_canonical(refine{i,j},mid,rot);
refine_j = map_points_to_canonical(refine{j,i},mid,rot);

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
if precomp
    Upf = -Uf_pair'*Npair;
else
    Upf = Uf_pair';
end
Ypf = Yf_pair;

DC = [];
YC = [];
C = [];
if opt.N_peanut
    rout_peanut = createPeanut(q_pair(1),q_pair(2),opt.N_peanut,0,R);
    rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
    proj_peanut = struct('project_charge',project_charge, ...
        'nsrc_c',[numel(rbase_in_c) numel(rbase_in_c)], ...
        'nsrc_f',[nsrc_f_i nsrc_f_j]);
    [DC,YC] = getPeanutBlockLaplace(rin_pair_c,rin_pair_f,rout_peanut,proj_peanut);

    if isfield(opt,'cmap') && opt.cmap
        C = -YC*(DC*Yf_pair*(Uf_pair'*Npair));
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
group.nout_f = numel(rout_base_f);
group.nsrc_f = [nsrc_f_i nsrc_f_j];
group.ntar_f = [ntar_i ntar_j];
group.rep_pair = [i j];
end

function z_canon = map_points_to_canonical(z,mid,rot)
z = z(:);
if isempty(z)
    z_canon = zeros(0,1);
    return
end

z_canon = conj(rot)*(z-mid);
end
