function [Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair_cache] = getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc)
%GETPAIRBASISSTOKES Build Stokes pair-basis pseudoinverse factors.
%
% Syntax:
%   [Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair_cache] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt)
%   [Uf,Yf,Up,Yp,Cmap,Cmap_FU,pair_cache] = getPairBasisStokes( ...
%       q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc)
%
% Anna Broms, Mar 2026

q = q(:);
P = opt.P;
N_peanut = getOptField(opt,'N_peanut',0);
use_pair_observable_map = logical(getOptField(opt,'cmap',false));
get_bndry_field = logical(getOptField(opt,'get_bndry_field',true));
reuse_pair_basis = logical(getOptField(opt,'reuse_pair_basis_by_sep',false));

shared_sep_tol = getOptField(opt,'shared_sep_tol',1e-2);


%Draw pair discretisation for debugging?
if isfield(opt,'pair_basis_debug') && ~isempty(opt.pair_basis_debug)
    debug = logical(opt.pair_basis_debug);
else
    debug = false;
end

%Show progress for building pair basis?
show_counter = opt.show_counter; 

% Need fine sources explicitly?
need_explicit_pair_sources = (N_peanut == 0) || ~use_pair_observable_map || get_bndry_field;

if reuse_pair_basis
    Uf = [];
    Yf = [];
    Up = [];
    Yp = [];
    Cmap = [];
    Cmap_FU = [];
else
    [Uf,Yf,Up,Yp,Cmap,Cmap_FU] = init_outputs(P,N_peanut,use_pair_observable_map,need_explicit_pair_sources);
end

pair_cache = init_pair_cache();
pair_cache.enabled = reuse_pair_basis;
pair_cache.shared_sep_tol = shared_sep_tol;

if isempty(pairs)
    return
end

N_f = opt.N_f;
pair_cache.meta = build_pair_meta(q,pairs,numel(rbase_in_c),N_f,opt);

if ~reuse_pair_basis
    total_pairs = size(pairs,1);
    processed_pairs = 0;

    for ii = 1:total_pairs
        i = pairs(ii,1);
        p2 = pairs(ii,2);

        pair = build_pair_data(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,ii,debug);
        if need_explicit_pair_sources
            Uf{i,p2} = pair.Upf;
            Yf{i,p2} = pair.Ypf;
        end
        if N_peanut && ~use_pair_observable_map
            Up{i,p2} = pair.DC;
            Yp{i,p2} = pair.YC;
        elseif N_peanut && use_pair_observable_map
            Cmap{i,p2} = pair.Cmap;
        end
        if use_pair_observable_map
            Cmap_FU{i,p2} = pair.Cmap_FU;
        end

        processed_pairs = processed_pairs + 1;
        if show_counter
            fprintf('getPairBasisStokes: processed pair %d/%d (%d,%d)\n', ...
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
total_pairs = size(pairs,1);
covered_pairs = 0;

for gg = 1:n_groups
    pair_cache.groups(gg) = build_pair_group(gg,rep_rows(gg),group_sep(gg), ...
        q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc, ...
        need_explicit_pair_sources,debug);

    if show_counter
        rep_i = pairs(rep_rows(gg),1);
        rep_j = pairs(rep_rows(gg),2);
        n_rep = sum(group_id == gg);
        covered_pairs = covered_pairs + n_rep;
        fprintf(['getPairBasisStokes: processed canonical group %d/%d ', ...
            'from pair (%d,%d), covers %d pairs -> %d/%d pairs covered\n'], ...
            gg,n_groups,rep_i,rep_j,n_rep,covered_pairs,total_pairs);
    end
end

for row = 1:size(pairs,1)
    gid = group_id(row);
    pair_cache.meta(row).group_id = gid;
    pair_cache.meta(row).sep = group_sep(gid);
end

if show_counter
    fprintf('getPairBasisStokes: reused %d canonical groups for %d pairs\n', ...
        n_groups,size(pairs,1));
end

end

function [Uf,Yf,Up,Yp,Cmap,Cmap_FU] = init_outputs(P,N_peanut,use_pair_observable_map,need_explicit_pair_sources)
if need_explicit_pair_sources
    Uf = cell(P);
    Yf = cell(P);
else
    Uf = [];
    Yf = [];
end

if N_peanut && ~use_pair_observable_map
    Up = cell(P);
    Yp = cell(P);
else
    Up = [];
    Yp = [];
end

if N_peanut && use_pair_observable_map
    Cmap = cell(P);
else
    Cmap = [];
end

if use_pair_observable_map
    Cmap_FU = cell(P);
else
    Cmap_FU = [];
end
end

function pair_cache = init_pair_cache()
pair_cache = struct();
pair_cache.enabled = false;
pair_cache.shared_sep_tol = [];
pair_cache.meta = repmat(struct('i',[],'j',[],'group_id',[],'sep',[], ...
    'mid',[],'rot',[],'phase_c',[],'phase_c_inv',[], ...
    'phase_f',[],'phase_f_inv',[]),0,1);
pair_cache.groups = repmat(init_pair_group(),0,1);
pair_cache.group_id = zeros(0,1);
pair_cache.group_sep = zeros(0,1);
pair_cache.representative_rows = zeros(0,1);
pair_cache.n_groups = 0;
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
    'phase_f',[],'phase_f_inv',[]),total_pairs,1);

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

function group = build_pair_group(group_id,row,group_sep,q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,need_explicit_pair_sources,debug)
pair = build_pair_data(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,row,debug,true);
group = init_pair_group();
group.group_id = group_id;
group.sep = group_sep;
group.q_pair = pair.q_pair;
group.rimage_canon = pair.rimage_canon;
group.refine_canon = pair.refine_canon;
group.Cmap = pair.Cmap;
group.Cmap_FU = pair.Cmap_FU;
group.Upair_colloc = pair.Upair_colloc;
group.Ucross_colloc = pair.Ucross_colloc;
group.Ecolloc = pair.Ecolloc;
group.rep_pair = pair.rep_pair;

if need_explicit_pair_sources
    group.Upf = pair.Upf;
    group.Ypf = pair.Ypf;
    group.DC = pair.DC;
    group.YC = pair.YC;
    group.Lf_pair = pair.Lf_pair;
    group.Lc_pair = pair.Lc_pair;
end
end

function pair = build_pair_data(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,row,debug,use_canonical)
if nargin < 11
    use_canonical = false;
end

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
pair.rin_pair = rin_pair;

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

[Uf_pair,Yf_pair] = getPairBlockStokes(rin_pair,rout_f,Lf_pair,Lr_pair);
Npair = evaluateCoarseOnPair(q_pair,rbase_in_c,rout_f);
Upf = -Uf_pair'*Npair; %

pair.Upf = Upf;
pair.Ypf = Yf_pair;
pair.Lf_pair = Lf_pair;
pair.Lc_pair = Lc_pair;

pair.DC = [];
pair.YC = [];
pair.Cmap = [];
pair.Cmap_FU = [];

if getOptField(opt,'N_peanut',0)
    rout_peanut = createPeanut(q_pair(1),q_pair(2),opt.N_peanut,0);
    rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
    [DC,YC] = getPeanutBlockStokes(rin_pair_c,rin_pair,rout_peanut,Lc_pair,Lf_pair);
    pair.DC = DC;
    pair.YC = YC;
    if logical(getOptField(opt,'cmap',false))
        pair.Cmap = -YC*(DC*Yf_pair*(Uf_pair'*Npair));
    end
else
    rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
end


if logical(getOptField(opt,'cmap',false))
    Kft_pair = getKftPair(Kf1,Kf2);
    pair.Cmap_FU = -Kft_pair*Yf_pair*(Uf_pair'*Npair);
end

pair.Upair_colloc = [];
pair.Ucross_colloc = [];
pair.Ecolloc = [];
if project_force
    rout_pair_c = [q_pair(1)+rout_base_c; q_pair(2)+rout_base_c];
    Ppair = eye(size(Lf_pair)) - Lf_pair;
    pair.Upair_colloc = stokSLPmat(rin_pair,rout_pair_c,1) * Ppair * Yf_pair * Upf;
    pair.Ecolloc = stokSLPmat(rin_pair_c,rout_pair_c,1);
    pair.Ucross_colloc = build_cross_pair_velocity_map(pair.Ecolloc,N_c,numel(rout_base_c));
end
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
