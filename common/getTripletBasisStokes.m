function triplet_cache = getTripletBasisStokes(q,rbase_in_c,rbase_in_f,rimage_pairs,refine,triplets,opt,svd_opts,pair_test_data)
%GETTRIPLETBASISSTOKES Build cached local 3-body Stokes correction bases.

if nargin < 8 || isempty(svd_opts)
    svd_opts = struct();
end
if nargin < 9
    pair_test_data = [];
end

q = q(:);
triplet_cache = init_triplet_cache();

if isempty(triplets)
    return
end

reuse_triplet_basis = logical(getOptField(opt,'reuse_triplet_basis_by_shape',true));
shape_tol = getOptField(opt,'shared_triplet_tol',1e-4);
N_c = opt.N_c;
N_f = opt.N_f;
debug = logical(getOptField(opt,'triplet_basis_debug',0));

triplet_cache.enabled = reuse_triplet_basis;
triplet_cache.shared_triplet_tol = shape_tol;
triplet_cache.meta = repmat(init_triplet_meta(),size(triplets,1),1);

for row = 1:size(triplets,1)
    triplet_cache.meta(row) = build_triplet_meta(q,triplets(row,:),N_c,N_f,opt);
end

group_desc = zeros(0,3);
group_rows = zeros(0,1);
group_id = zeros(size(triplets,1),1);

for row = 1:size(triplets,1)
    desc = triplet_cache.meta(row).descriptor(:).';
    if reuse_triplet_basis
        gid = find_group(group_desc,desc,shape_tol);
    else
        gid = 0;
    end

    if gid == 0
        gid = size(group_desc,1) + 1;
        group_desc(gid,:) = desc;
        group_rows(gid,1) = row;
    end
    group_id(row) = gid;
    triplet_cache.meta(row).group_id = gid;
end

triplet_cache.group_descriptor = group_desc;
triplet_cache.representative_rows = group_rows;
triplet_cache.n_groups = size(group_desc,1);
triplet_cache.groups = repmat(init_triplet_group(),triplet_cache.n_groups,1);

for gid = 1:triplet_cache.n_groups
    rep_row = group_rows(gid);
    triplet_cache.groups(gid) = build_triplet_group( ...
        gid,triplet_cache.meta(rep_row),q,rbase_in_c,rbase_in_f, ...
        rimage_pairs,refine,opt,svd_opts,pair_test_data);
end

if debug
    fprintf('getTripletBasisStokes: running triplet LSQ debug for %d group(s)\n', ...
        triplet_cache.n_groups);
end

end

function cache = init_triplet_cache()
cache = struct();
cache.enabled = false;
cache.shared_triplet_tol = [];
cache.meta = repmat(init_triplet_meta(),0,1);
cache.groups = repmat(init_triplet_group(),0,1);
cache.group_descriptor = zeros(0,3);
cache.representative_rows = zeros(0,1);
cache.n_groups = 0;
end

function meta = init_triplet_meta()
meta = struct('ids_sorted',[],'ids',[],'group_id',[],'descriptor',[], ...
    'anchor',[],'rot',[],'q_canon',[],'phase_c',[],'phase_c_inv',[], ...
    'phase_f',[],'phase_f_inv',[]);
end

function group = init_triplet_group()
group = struct('group_id',[],'descriptor',[],'q_trip',[],'rimage_canon',{{}}, ...
    'refine_canon',{{}},'Upf',[],'Ypf',[],'Upf_full',[],'Ypf_full',[], ...
    'lsq_mode','','rep_ids',[]);
end

function meta = build_triplet_meta(q,triplet,nc,nf,opt)
ids_sorted = sort(triplet(:)).';
z = q(ids_sorted);
perm_list = perms(1:3);

best_key = [];
best_perm = [];
best_q_canon = [];
best_rot = [];
best_anchor = [];

for kk = 1:size(perm_list,1)
    perm_try = perm_list(kk,:);
    zloc = z(perm_try);
    edge12 = zloc(2) - zloc(1);
    d12 = abs(edge12);
    if d12 == 0
        continue
    end

    rot = edge12/d12;
    z3 = conj(rot)*(zloc(3) - zloc(1));
    if imag(z3) < -1e-12
        continue
    end

    key = [d12, real(z3), imag(z3), perm_try];
    if isempty(best_key) || lexicographically_smaller(key,best_key)
        best_key = key;
        best_perm = perm_try;
        best_q_canon = [0; d12; z3];
        best_rot = rot;
        best_anchor = zloc(1);
    end
end

if isempty(best_perm)
    perm_try = [1 2 3];
    zloc = z(perm_try);
    edge12 = zloc(2) - zloc(1);
    d12 = abs(edge12);
    rot = edge12/max(d12,eps);
    z3 = conj(rot)*(zloc(3) - zloc(1));
    best_perm = perm_try;
    best_q_canon = [0; d12; z3];
    best_rot = rot;
    best_anchor = zloc(1);
end

meta = init_triplet_meta();
meta.ids_sorted = ids_sorted;
meta.ids = ids_sorted(best_perm);
meta.anchor = best_anchor;
meta.rot = best_rot;
meta.q_canon = best_q_canon;
meta.descriptor = [best_q_canon(2), real(best_q_canon(3)), imag(best_q_canon(3))];
meta.phase_c = getUniformCircleRotationSpec(nc,best_rot,opt);
meta.phase_c_inv = invertUniformCircleRotationSpec(meta.phase_c);
meta.phase_f = getUniformCircleRotationSpec(nf,best_rot,opt);
meta.phase_f_inv = invertUniformCircleRotationSpec(meta.phase_f);
end

function tf = lexicographically_smaller(a,b)
tf = false;
for k = 1:min(numel(a),numel(b))
    if a(k) < b(k)
        tf = true;
        return
    elseif a(k) > b(k)
        return
    end
end
tf = numel(a) < numel(b);
end

function gid = find_group(group_desc,desc,tol)
gid = 0;
for kk = 1:size(group_desc,1)
    if max(abs(group_desc(kk,:) - desc)) <= tol
        gid = kk;
        return
    end
end
end

function group = build_triplet_group(gid,meta,q,rbase_in_c,rbase_in_f,rimage_pairs,refine,opt,svd_opts,pair_test_data)
group = init_triplet_group();
group.group_id = gid;
group.descriptor = meta.descriptor;
group.q_trip = meta.q_canon;
group.rep_ids = meta.ids;

[rimage_canon,refine_canon] = collect_triplet_geometry(meta,q,rimage_pairs,refine);
group.rimage_canon = rimage_canon;
group.refine_canon = refine_canon;

[rin_trip,rout_trip,rout_cells] = build_triplet_point_sets(group.q_trip,rbase_in_f,refine_canon,rimage_canon,opt);
K1 = [];
K2 = [];
K3 = [];

if logical(getOptField(opt,'project_force',false))
    K1 = getKmat2D([group.q_trip(1)+rbase_in_f; rimage_canon{1}],group.q_trip(1));
    K2 = getKmat2D([group.q_trip(2)+rbase_in_f; rimage_canon{2}],group.q_trip(2));
    K3 = getKmat2D([group.q_trip(3)+rbase_in_f; rimage_canon{3}],group.q_trip(3));

    B1 = getKmat2D(rout_cells{1},group.q_trip(1));
    B2 = getKmat2D(rout_cells{2},group.q_trip(2));
    B3 = getKmat2D(rout_cells{3},group.q_trip(3));

    Lf_trip = getLfTriplet(K1,K2,K3);
    Lr_trip = getLrTriplet(B1,B2,B3,K1,K2,K3);
else
    Lf_trip = [];
    Lr_trip = [];
end

[Uf_trip,Yf_trip,Atot] = getTripletBlockStokes(rin_trip,rout_trip,Lf_trip,Lr_trip,svd_opts);
Ntrip = evaluateCoarseOnTriplet(group.q_trip,rbase_in_c,rout_cells);
group.Upf_full = -Uf_trip'*Ntrip;
group.Ypf_full = Yf_trip;
group.lsq_mode = get_triplet_solver_mode(opt);

pair_debug = [];
need_pair_maps = triplet_lsq_mode_needs_pair_maps(group.lsq_mode) || ...
    (logical(getOptField(opt,'triplet_basis_debug',0)) && ...
    triplet_debug_needs_pair_maps(opt));
if need_pair_maps
    pair_debug = build_triplet_pair_debug_maps(meta.ids,group.q_trip, ...
        cellfun(@numel,group.rimage_canon(:)).',opt,pair_test_data);
end

[group.Upf,group.Ypf] = build_triplet_solver_factors(group.lsq_mode, ...
    Uf_trip,Yf_trip,Atot,Ntrip,pair_debug,opt,svd_opts,meta.ids);

if logical(getOptField(opt,'triplet_basis_debug',0))
    run_triplet_lsq_debug_test(meta.ids,gid,group.q_trip,rbase_in_c,rin_trip, ...
        Uf_trip,Yf_trip,Atot,Lf_trip,K1,K2,K3,rout_cells,Ntrip, ...
        cellfun(@numel,group.rimage_canon(:)).',opt,svd_opts,pair_debug);
end
end

function [rimage_canon,refine_canon] = collect_triplet_geometry(meta,q,rimage_pairs,refine)
rimage_canon = cell(3,1);
refine_canon = cell(3,1);
ids = meta.ids;

for body = 1:3
    neigh = setdiff(1:3,body,'stable');
    id_body = ids(body);
    id_nb1 = ids(neigh(1));
    id_nb2 = ids(neigh(2));

    rimage_actual = [rimage_pairs{id_body,id_nb1}; rimage_pairs{id_body,id_nb2}];
    refine_actual = [refine{id_body,id_nb1}; refine{id_body,id_nb2}];

    rimage_canon{body} = map_points_to_canonical(rimage_actual,meta);
    refine_canon{body} = map_points_to_canonical(refine_actual,meta);
end
end

function z_canon = map_points_to_canonical(z_actual,meta)
z_actual = z_actual(:);
if isempty(z_actual)
    z_canon = zeros(0,1);
    return
end

z_canon = conj(meta.rot)*(z_actual - meta.anchor);
end

function [rin_trip,rout_trip,rout_cells] = build_triplet_point_sets(q_trip,rbase_in_f,refine_cells,rimage_cells,opt)
nout = ceil(opt.a_f*opt.N_f);
t = linspace(0,2*pi,nout+1)';
t = t(1:end-1);
rout_base = cos(t) + 1i*sin(t);

rin_cells = cell(3,1);
rout_cells = cell(3,1);
for body = 1:3
    rin_cells{body} = [q_trip(body) + rbase_in_f; rimage_cells{body}];
    rout_cells{body} = [q_trip(body) + rout_base; refine_cells{body}];
end

rin_trip = vertcat(rin_cells{:});
rout_trip = vertcat(rout_cells{:});
end

function run_triplet_lsq_debug_test(ids,gid,q_trip,rbase_in_c,rin_trip,Uf_trip,Yf_trip,Atot, ...
    Lf_trip,K1,K2,K3,rout_cells,Ntrip_colloc,im_counts_trip,opt,svd_opts,pair_debug)
N_c = opt.N_c;
[lambda_1,lambda_2,lambda_3] = build_triplet_debug_rhs(N_c,gid,opt);
rhs_trip = [lambda_1(1:N_c); lambda_2(1:N_c); lambda_3(1:N_c); ...
    lambda_1(N_c+1:end); lambda_2(N_c+1:end); lambda_3(N_c+1:end)];
mode_list = get_triplet_debug_mode_list(opt);

rout_test_cells = build_gap_refined_triplet_boundary_nodes(q_trip,2*cellfun(@numel,rout_cells(:)));
rout_test = vertcat(rout_test_cells{:});
S_test = stokSLPmat(rin_trip,rout_test,1);
if isempty(Lf_trip)
    A_test = S_test;
else
    B1_test = getKmat2D(rout_test_cells{1},q_trip(1));
    B2_test = getKmat2D(rout_test_cells{2},q_trip(2));
    B3_test = getKmat2D(rout_test_cells{3},q_trip(3));
    Lr_test = getLrTriplet(B1_test,B2_test,B3_test,K1,K2,K3);
    A_test = S_test - S_test*Lf_trip + Lr_test;
end

Ntrip_test = evaluateCoarseOnTriplet(q_trip,rbase_in_c,rout_test_cells);
r_trip_test_full = -Ntrip_test*rhs_trip;

r_trip_colloc_full = -Ntrip_colloc*rhs_trip;
beta_full = Yf_trip*(Uf_trip'*r_trip_colloc_full);

need_pair_debug = any(strcmp(mode_list,'rhs_subtracted')) || ...
    any(strcmp(mode_list,'pair_complement')) || ...
    any(strcmp(mode_list,'rhs_subtracted_pair_complement'));
if need_pair_debug
    if isempty(pair_debug) || ~isstruct(pair_debug) || ~pair_debug.available
        fprintf(['getTripletBasisStokes debug triplet (%d,%d,%d): ', ...
            'pair-based LSQ experiments requested, but pair debug data is unavailable. ', ...
            'Running full mode only.\n'],ids(1),ids(2),ids(3));
        mode_list = {'full'};
    else
        pair_debug.beta_pair_sum = pair_debug.Bsum*rhs_trip;
        pair_debug.r_rhs_colloc = r_trip_colloc_full - Atot*pair_debug.beta_pair_sum;
        pair_debug.r_rhs_test = r_trip_test_full - A_test*pair_debug.beta_pair_sum;
        pair_debug.pair_removed_abs = norm(Atot*pair_debug.beta_pair_sum + Ntrip_colloc*rhs_trip);
        pair_debug.pair_removed_rel = pair_debug.pair_removed_abs / ...
            max(norm(Ntrip_colloc*rhs_trip),eps);
        ABpair = Atot*pair_debug.Bpair;
        pair_debug.Q = orth(ABpair);
        if isempty(pair_debug.Q)
            pair_debug.r_pair_frac = 0;
            pair_debug.r_comp_frac = 1;
        else
            r_pair_proj = pair_debug.Q*(pair_debug.Q'*r_trip_colloc_full);
            pair_debug.r_pair_frac = norm(r_pair_proj)/max(norm(r_trip_colloc_full),eps);
            pair_debug.r_comp_frac = norm(r_trip_colloc_full-r_pair_proj)/ ...
                max(norm(r_trip_colloc_full),eps);
        end
        fprintf(['getTripletBasisStokes debug triplet (%d,%d,%d): ', ...
            '||A*beta_pair_sum + Ntrip*rhs|| = %.3e (rel %.3e), rank(Bpair) = %d, dim(Z) = %d\n'], ...
            ids(1),ids(2),ids(3),pair_debug.pair_removed_abs,pair_debug.pair_removed_rel, ...
            pair_debug.rank_bpair,pair_debug.dim_z);
        fprintf(['getTripletBasisStokes debug triplet (%d,%d,%d): ', ...
            '||QQ^T r||/||r|| = %.3e, ||(I-QQ^T)r||/||r|| = %.3e\n'], ...
            ids(1),ids(2),ids(3),pair_debug.r_pair_frac,pair_debug.r_comp_frac);
        if pair_debug.dim_z == 0
            fprintf(['getTripletBasisStokes debug triplet (%d,%d,%d): ', ...
                'WARNING: Z is empty, so pair-complement tests will return beta = 0.\n'], ...
                ids(1),ids(2),ids(3));
        end
    end
end

for mode_it = 1:numel(mode_list)
    mode_name = mode_list{mode_it};
    [beta_mode,r_colloc,r_test,leakage] = solve_triplet_debug_mode( ...
        mode_name,beta_full,r_trip_colloc_full,r_trip_test_full,Atot,A_test,Uf_trip,Yf_trip, ...
        pair_debug,svd_opts);
    [abs_off,abs_pw_off,rel_off,rel_pw_off,abs_coll,abs_pw_coll,rel_coll,rel_pw_coll] = ...
        evaluate_triplet_debug_solution(A_test,Atot,beta_mode,r_test,r_colloc);

    fig_id = 9100 + 8*(gid-1) + mode_it;
    plot_abs_rel_error_compare(fig_id,abs_pw_off,abs_pw_coll,rel_pw_off,rel_pw_coll, ...
        'Boundary node number',sprintf(['Triplet (%d,%d,%d) %s LSQ: ', ...
        'abs_{inf} off/coll = %.3e / %.3e, rel_{inf} off/coll = %.3e / %.3e'], ...
        ids(1),ids(2),ids(3),mode_name,abs_off,abs_coll,rel_off,rel_coll));

    fprintf(['getTripletBasisStokes debug triplet (%d,%d,%d) mode=%s: ', ...
        'absinf off/coll = %.3e / %.3e, relinf off/coll = %.3e / %.3e, ||beta||_2 = %.3e'], ...
        ids(1),ids(2),ids(3),mode_name,abs_off,abs_coll,rel_off,rel_coll,norm(beta_mode));
    if ~isnan(leakage)
        fprintf(', leakage = %.3e',leakage);
    end
    fprintf('\n');
end
end

function mode_name = get_triplet_solver_mode(opt)
mode_name = lower(char(getOptField(opt,'triplet_lsq_mode','full')));
valid_modes = {'full','rhs_subtracted','pair_complement','rhs_subtracted_pair_complement'};
if ~any(strcmp(mode_name,valid_modes))
    fprintf('getTripletBasisStokes: unknown triplet_lsq_mode "%s", using full.\n', ...
        mode_name);
    mode_name = 'full';
end
end

function tf = triplet_lsq_mode_needs_pair_maps(mode_name)
tf = any(strcmp(mode_name,{'rhs_subtracted','pair_complement', ...
    'rhs_subtracted_pair_complement'}));
end

function tf = triplet_debug_needs_pair_maps(opt)
mode_list = get_triplet_debug_mode_list(opt);
tf = any(strcmp(mode_list,'rhs_subtracted')) || any(strcmp(mode_list,'pair_complement')) || ...
    any(strcmp(mode_list,'rhs_subtracted_pair_complement'));
end

function [Upf_sel,Ypf_sel] = build_triplet_solver_factors(mode_name,Uf_trip,Yf_trip,Atot,Ntrip, ...
    pair_debug,opt,svd_opts,ids)
Upf_full = -Uf_trip'*Ntrip;
Ypf_full = Yf_trip;

switch mode_name
    case 'full'
        Upf_sel = Upf_full;
        Ypf_sel = Ypf_full;

    case 'rhs_subtracted'
        require_pair_maps_for_triplet_mode(pair_debug,mode_name,ids);
        rhs_map = -Ntrip - Atot*pair_debug.Bsum;
        Upf_sel = Uf_trip'*rhs_map;
        Ypf_sel = Ypf_full;

    case 'pair_complement'
        require_pair_maps_for_triplet_mode(pair_debug,mode_name,ids);
        [Upf_sel,Ypf_sel] = build_triplet_complement_factors( ...
            Atot,-Ntrip,pair_debug,opt,svd_opts,ids,mode_name);

    case 'rhs_subtracted_pair_complement'
        require_pair_maps_for_triplet_mode(pair_debug,mode_name,ids);
        rhs_map = -Ntrip - Atot*pair_debug.Bsum;
        [Upf_sel,Ypf_sel] = build_triplet_complement_factors( ...
            Atot,rhs_map,pair_debug,opt,svd_opts,ids,mode_name);

    otherwise
        error('Unsupported triplet_lsq_mode "%s".',mode_name);
end
end

function require_pair_maps_for_triplet_mode(pair_debug,mode_name,ids)
if isempty(pair_debug) || ~isstruct(pair_debug) || ~pair_debug.available
    error(['getTripletBasisStokes: triplet_lsq_mode "%s" for triplet (%d,%d,%d) ', ...
        'requires pair basis data, but none was provided.'], ...
        mode_name,ids(1),ids(2),ids(3));
end
end

function [Upf_sel,Ypf_sel] = build_triplet_complement_factors(Atot,rhs_map,pair_debug,opt,svd_opts,ids,mode_name)
Z = pair_debug.Z;
if isempty(Z)
    fprintf(['getTripletBasisStokes: triplet_lsq_mode "%s" for triplet (%d,%d,%d) ', ...
        'has empty complement basis Z, returning beta = 0.\n'], ...
        mode_name,ids(1),ids(2),ids(3));
    Upf_sel = zeros(0,size(rhs_map,2));
    Ypf_sel = zeros(size(Atot,2),0);
    return
end

r_eff = project_rhs_off_pair_field(Atot,pair_debug.Bpair,rhs_map);
[Yz,Uz] = getPseudoFactors(Atot*Z,1e-9,0,svd_opts);
Upf_sel = Uz'*r_eff;
Ypf_sel = Z*Yz;
end

function mode_list = get_triplet_debug_mode_list(opt)
mode_name = lower(char(getOptField(opt,'triplet_lsq_test_mode','full')));
switch mode_name
    case 'full'
        mode_list = {'full'};
    case 'rhs_subtracted'
        mode_list = {'rhs_subtracted'};
    case 'pair_complement'
        mode_list = {'pair_complement'};
    case 'rhs_subtracted_pair_complement'
        mode_list = {'rhs_subtracted_pair_complement'};
    case 'compare_all'
        mode_list = {'full','rhs_subtracted','pair_complement', ...
            'rhs_subtracted_pair_complement'};
    otherwise
        fprintf('getTripletBasisStokes: unknown triplet_lsq_test_mode "%s", using full.\n', ...
            mode_name);
        mode_list = {'full'};
end
end

function pair_debug = build_triplet_pair_debug_maps(ids,q_trip,im_counts_expected,opt,pair_test_data)
pair_debug = struct();
pair_debug.available = false;
pair_debug.B12 = [];
pair_debug.B13 = [];
pair_debug.B23 = [];
pair_debug.Bsum = [];
pair_debug.Bpair = [];
pair_debug.Z = [];
pair_debug.Q = [];
pair_debug.rank_bpair = 0;
pair_debug.dim_z = 0;
pair_debug.beta_pair_sum = [];
pair_debug.r_rhs_colloc = [];
pair_debug.r_rhs_test = [];
pair_debug.pair_removed_abs = NaN;
pair_debug.pair_removed_rel = NaN;
pair_debug.r_pair_frac = NaN;
pair_debug.r_comp_frac = NaN;

if isempty(pair_test_data) || ~isstruct(pair_test_data) || ...
        ~isfield(pair_test_data,'pair_lookup') || ~isfield(pair_test_data,'pair_cache')
    return
end

pair_lookup = pair_test_data.pair_lookup;
pair_cache = pair_test_data.pair_cache;
if isempty(pair_cache) || ~isfield(pair_cache,'groups') || isempty(pair_cache.groups)
    return
end

edge_list = [1 2; 1 3; 2 3];
directed_counts = zeros(3,3);
edge_data = repmat(struct('a',[],'b',[],'pair_row',[],'group',[],'edge_rot',[], ...
    'phase_c',[],'phase_f_inv',[],'im_a',[],'im_b',[]),size(edge_list,1),1);
for edge = 1:size(edge_list,1)
    local_a = edge_list(edge,1);
    local_b = edge_list(edge,2);
    pair_row = pair_lookup(ids(local_a),ids(local_b));
    if pair_row == 0
        return
    end
    meta_pair = pair_cache.meta(pair_row);
    group_pair = pair_cache.groups(meta_pair.group_id);
    edge_vec = q_trip(local_b) - q_trip(local_a);
    edge_sep = abs(edge_vec);
    if edge_sep == 0
        edge_rot = 1;
    else
        edge_rot = edge_vec/edge_sep;
    end

    im_a = numel(group_pair.rimage_canon{1});
    im_b = numel(group_pair.rimage_canon{2});
    directed_counts(local_a,local_b) = im_a;
    directed_counts(local_b,local_a) = im_b;

    edge_data(edge).a = local_a;
    edge_data(edge).b = local_b;
    edge_data(edge).pair_row = pair_row;
    edge_data(edge).group = group_pair;
    edge_data(edge).edge_rot = edge_rot;
    edge_data(edge).phase_c = getUniformCircleRotationSpec(opt.N_c,edge_rot,opt);
    phase_f = getUniformCircleRotationSpec(opt.N_f,edge_rot,opt);
    edge_data(edge).phase_f_inv = invertUniformCircleRotationSpec(phase_f);
    edge_data(edge).im_a = im_a;
    edge_data(edge).im_b = im_b;
end

im_counts_trip = sum(directed_counts,2).';
if any(im_counts_trip ~= im_counts_expected)
    fprintf(['getTripletBasisStokes debug triplet (%d,%d,%d): ', ...
        'WARNING: reconstructed triplet image counts [%s] do not match cached counts [%s].\n'], ...
        ids(1),ids(2),ids(3),num2str(im_counts_trip),num2str(im_counts_expected));
end
idx_trip = get_triplet_source_indices_debug(opt.N_f,im_counts_trip);

B_edge = cell(3,1);
M_edge = cell(3,1);
for edge = 1:size(edge_list,1)
    M_edge{edge} = build_triplet_edge_source_map(edge_data(edge),idx_trip,directed_counts,opt);
    P_edge = build_triplet_edge_selector_map(edge_data(edge),opt.N_c);
    B_edge{edge} = M_edge{edge}*P_edge;
end

pair_debug.B12 = B_edge{1};
pair_debug.B13 = B_edge{2};
pair_debug.B23 = B_edge{3};
pair_debug.Bsum = B_edge{1} + B_edge{2} + B_edge{3};
pair_debug.Bpair = [M_edge{1}, M_edge{2}, M_edge{3}];
pair_debug.rank_bpair = rank(pair_debug.Bpair,opt.triplet_pair_subspace_tol);
pair_debug.Z = null(pair_debug.Bpair.',opt.triplet_pair_subspace_tol);
pair_debug.dim_z = size(pair_debug.Z,2);
pair_debug.available = true;
end

function M_edge = build_triplet_edge_source_map(edge_data,idx_trip,directed_counts,opt)
pair_map = edge_data.group.Ypf*edge_data.group.Upf;
M_edge = zeros(idx_trip.ntot,size(pair_map,2));
for col = 1:size(pair_map,2)
    beta_pair = rotateStokesPairSourceVector(pair_map(:,col),opt.N_f, ...
        edge_data.im_a,edge_data.im_b,edge_data.phase_f_inv,edge_data.edge_rot);
    M_edge(:,col) = embed_pair_vector_into_triplet_debug(beta_pair,edge_data, ...
        idx_trip,directed_counts,opt.N_f);
end
end

function P_edge = build_triplet_edge_selector_map(edge_data,N_c)
P_edge = zeros(4*N_c,6*N_c);
for col = 1:6*N_c
    rhs_trip = zeros(6*N_c,1);
    rhs_trip(col) = 1;
    rhs_pair = select_triplet_pair_coarse(rhs_trip,edge_data.a,edge_data.b,N_c);
    P_edge(:,col) = rotatePairOrderedStokesData(rhs_pair,N_c, ...
        edge_data.phase_c,conj(edge_data.edge_rot));
end
end

function beta_embed = embed_pair_vector_into_triplet_debug(beta_pair,edge_data,idx_trip,directed_counts,N_f)
pair_idx = get_pair_source_indices_debug(N_f,edge_data.im_a,edge_data.im_b);
beta_embed = zeros(idx_trip.ntot,1);

beta_embed(idx_trip_body_fine(idx_trip,edge_data.a,'x')) = beta_pair(pair_idx.f1x);
beta_embed(idx_trip_body_fine(idx_trip,edge_data.a,'y')) = beta_pair(pair_idx.f1y);
beta_embed(idx_trip_body_fine(idx_trip,edge_data.b,'x')) = beta_pair(pair_idx.f2x);
beta_embed(idx_trip_body_fine(idx_trip,edge_data.b,'y')) = beta_pair(pair_idx.f2y);

[extra_idx_ax,extra_idx_ay] = get_triplet_edge_chunk_indices_debug(idx_trip,directed_counts, ...
    edge_data.a,edge_data.b);
[extra_idx_bx,extra_idx_by] = get_triplet_edge_chunk_indices_debug(idx_trip,directed_counts, ...
    edge_data.b,edge_data.a);

beta_embed(extra_idx_ax) = beta_pair(pair_idx.e1x);
beta_embed(extra_idx_ay) = beta_pair(pair_idx.e1y);
beta_embed(extra_idx_bx) = beta_pair(pair_idx.e2x);
beta_embed(extra_idx_by) = beta_pair(pair_idx.e2y);
end

function rhs_pair = select_triplet_pair_coarse(rhs_trip,local_a,local_b,N_c)
rhs_pair = zeros(4*N_c,1);
rhs_pair(1:N_c) = rhs_trip((local_a-1)*N_c+1:local_a*N_c);
rhs_pair(N_c+1:2*N_c) = rhs_trip((local_b-1)*N_c+1:local_b*N_c);
rhs_pair(2*N_c+1:3*N_c) = rhs_trip(3*N_c + (local_a-1)*N_c + 1:3*N_c + local_a*N_c);
rhs_pair(3*N_c+1:4*N_c) = rhs_trip(3*N_c + (local_b-1)*N_c + 1:3*N_c + local_b*N_c);
end

function [beta_mode,r_colloc,r_test,leakage] = solve_triplet_debug_mode( ...
    mode_name,beta_full,r_trip_colloc_full,r_trip_test_full,A_colloc,A_test,Uf_trip,Yf_trip, ...
    pair_debug,svd_opts)
leakage = NaN;
switch mode_name
    case 'full'
        beta_mode = beta_full;
        r_colloc = r_trip_colloc_full;
        r_test = r_trip_test_full;

    case 'rhs_subtracted'
        beta_mode = Yf_trip*(Uf_trip'*pair_debug.r_rhs_colloc);
        r_colloc = pair_debug.r_rhs_colloc;
        r_test = pair_debug.r_rhs_test;

    case 'pair_complement'
        r_colloc = project_rhs_off_pair_field(A_colloc,pair_debug.Bpair,r_trip_colloc_full);
        r_test = project_rhs_off_pair_field(A_test,pair_debug.Bpair,r_trip_test_full);
        beta_mode = solve_triplet_complement_lsq(A_colloc,pair_debug.Z,r_colloc,svd_opts);
        leakage = pair_span_leakage(pair_debug.Bpair,beta_mode);

    case 'rhs_subtracted_pair_complement'
        r_colloc = project_rhs_off_pair_field(A_colloc,pair_debug.Bpair,pair_debug.r_rhs_colloc);
        r_test = project_rhs_off_pair_field(A_test,pair_debug.Bpair,pair_debug.r_rhs_test);
        beta_mode = solve_triplet_complement_lsq(A_colloc,pair_debug.Z,r_colloc,svd_opts);
        leakage = pair_span_leakage(pair_debug.Bpair,beta_mode);

    otherwise
        error('Unknown triplet debug mode "%s".',mode_name);
end
end

function rhs_proj = project_rhs_off_pair_field(Amap,Bpair,rhs_vec)
if isempty(rhs_vec) || isempty(Bpair)
    rhs_proj = rhs_vec;
    return
end

Qpair = orth(Amap*Bpair);
if isempty(Qpair)
    rhs_proj = rhs_vec;
    return
end

rhs_proj = rhs_vec - Qpair*(Qpair'*rhs_vec);
end

function beta = solve_triplet_complement_lsq(A_colloc,Z,rhs_vec,svd_opts)
if isempty(Z)
    beta = zeros(size(A_colloc,2),1);
    return
end

[Yz,Uz] = getPseudoFactors(A_colloc*Z,1e-9,0,svd_opts);
gamma = Yz*(Uz'*rhs_vec);
beta = Z*gamma;
end

function leakage = pair_span_leakage(Bpair,beta)
if isempty(beta)
    leakage = 0;
    return
end
leakage = norm(Bpair.'*beta)/max(norm(beta),eps);
end

function [abs_off,abs_pw_off,rel_off,rel_pw_off,abs_coll,abs_pw_coll,rel_coll,rel_pw_coll] = ...
    evaluate_triplet_debug_solution(A_test,A_colloc,beta,r_test,r_colloc)
res_test = A_test*beta - r_test;
[abs_off,abs_pw_off,rel_off,rel_pw_off] = triplet_error_metrics(res_test,r_test);

res_coll = A_colloc*beta - r_colloc;
[abs_coll,abs_pw_coll,rel_coll,rel_pw_coll] = triplet_error_metrics(res_coll,r_colloc);
end

function idx = get_triplet_source_indices_debug(Nf,im_counts)
im1 = im_counts(1);
im2 = im_counts(2);
im3 = im_counts(3);

idx = struct();
idx.f1x = 1:Nf;
idx.e1x = Nf+1:Nf+im1;
idx.f2x = Nf+im1+1:2*Nf+im1;
idx.e2x = 2*Nf+im1+1:2*Nf+im1+im2;
idx.f3x = 2*Nf+im1+im2+1:3*Nf+im1+im2;
idx.e3x = 3*Nf+im1+im2+1:3*Nf+im1+im2+im3;

base_y = 3*Nf + im1 + im2 + im3;
idx.f1y = base_y + (1:Nf);
idx.e1y = base_y + Nf + (1:im1);
idx.f2y = base_y + Nf + im1 + (1:Nf);
idx.e2y = base_y + 2*Nf + im1 + (1:im2);
idx.f3y = base_y + 2*Nf + im1 + im2 + (1:Nf);
idx.e3y = base_y + 3*Nf + im1 + im2 + (1:im3);
idx.ntot = 6*Nf + 2*(im1 + im2 + im3);
end

function idx = get_pair_source_indices_debug(Nf,im_a,im_b)
idx = struct();
idx.f1x = 1:Nf;
idx.e1x = Nf+1:Nf+im_a;
idx.f2x = Nf+im_a+1:2*Nf+im_a;
idx.e2x = 2*Nf+im_a+1:2*Nf+im_a+im_b;
idx.f1y = 2*Nf+im_a+im_b+1:3*Nf+im_a+im_b;
idx.e1y = 3*Nf+im_a+im_b+1:3*Nf+2*im_a+im_b;
idx.f2y = 3*Nf+2*im_a+im_b+1:4*Nf+2*im_a+im_b;
idx.e2y = 4*Nf+2*im_a+im_b+1:4*Nf+2*im_a+2*im_b;
end

function idx_body = idx_trip_body_fine(idx_trip,body,coord)
switch body
    case 1
        idx_x = idx_trip.f1x;
        idx_y = idx_trip.f1y;
    case 2
        idx_x = idx_trip.f2x;
        idx_y = idx_trip.f2y;
    case 3
        idx_x = idx_trip.f3x;
        idx_y = idx_trip.f3y;
    otherwise
        error('Invalid triplet body index.');
end

if coord == 'x'
    idx_body = idx_x;
else
    idx_body = idx_y;
end
end

function [extra_idx_x,extra_idx_y] = get_triplet_edge_chunk_indices_debug(idx_trip,directed_counts,body,neigh)
others = setdiff(1:3,body,'stable');
if neigh == others(1)
    offset = 0;
    count = directed_counts(body,others(1));
elseif neigh == others(2)
    offset = directed_counts(body,others(1));
    count = directed_counts(body,others(2));
else
    error('Invalid triplet neighbour index.');
end

switch body
    case 1
        body_x = idx_trip.e1x;
        body_y = idx_trip.e1y;
    case 2
        body_x = idx_trip.e2x;
        body_y = idx_trip.e2y;
    case 3
        body_x = idx_trip.e3x;
        body_y = idx_trip.e3y;
    otherwise
        error('Invalid triplet body index.');
end

if count == 0
    extra_idx_x = zeros(0,1);
    extra_idx_y = zeros(0,1);
else
    extra_idx_x = body_x(offset+1:offset+count);
    extra_idx_y = body_y(offset+1:offset+count);
end
end

function [lambda_1,lambda_2,lambda_3] = build_triplet_debug_rhs(N_c,row,opt)
has_lam1 = isfield(opt,'triplet_basis_lambda_1') && ~isempty(opt.triplet_basis_lambda_1);
has_lam2 = isfield(opt,'triplet_basis_lambda_2') && ~isempty(opt.triplet_basis_lambda_2);
has_lam3 = isfield(opt,'triplet_basis_lambda_3') && ~isempty(opt.triplet_basis_lambda_3);
if has_lam1 || has_lam2 || has_lam3
    if ~(has_lam1 && has_lam2 && has_lam3)
        error(['triplet_basis_lambda_1, triplet_basis_lambda_2, and ', ...
            'triplet_basis_lambda_3 must be set together when debug is enabled.']);
    end
    lambda_1 = opt.triplet_basis_lambda_1(:);
    lambda_2 = opt.triplet_basis_lambda_2(:);
    lambda_3 = opt.triplet_basis_lambda_3(:);
    if numel(lambda_1) ~= 2*N_c || numel(lambda_2) ~= 2*N_c || numel(lambda_3) ~= 2*N_c
        error('Each triplet_basis_lambda_* vector must have length 2*N_c.');
    end
    return
end

t = linspace(0,2*pi,N_c+1)';
t = t(1:end-1);
phi = 0.29*row;

lambda1_x = cos(t+phi) + 0.25*sin(2*t-0.2) + 0.18*cos(3*t+0.4);
lambda1_y = 0.7*lambda1_x - 0.15*sin(3*t+phi);
lambda2_x = sin(t-0.3*phi) - 0.2*cos(2*t+0.7) + 0.11*sin(4*t-0.9);
lambda2_y = -0.55*lambda2_x + 0.17*cos(2*t-phi);
lambda3_x = cos(1.4*t+0.2*phi) + 0.22*sin(3*t-0.5) - 0.13*cos(5*t+0.1);
lambda3_y = -0.4*lambda3_x + 0.19*sin(2*t+0.6*phi);

lambda_1 = [lambda1_x; lambda1_y];
lambda_2 = [lambda2_x; lambda2_y];
lambda_3 = [lambda3_x; lambda3_y];
end

function rout_cells = build_gap_refined_triplet_boundary_nodes(q_trip,nout_vec)
rout_cells = cell(3,1);
for body = 1:3
    nout = max(16,round(nout_vec(body)));
    others = setdiff(1:3,body,'stable');
    theta_1 = angle(q_trip(others(1)) - q_trip(body));
    theta_2 = angle(q_trip(others(2)) - q_trip(body));

    n_uniform = max(12,round(0.45*nout));
    n_uniform = min(n_uniform,nout);
    n_cluster = nout - n_uniform;
    n_cluster_1 = floor(n_cluster/2);
    n_cluster_2 = n_cluster - n_cluster_1;

    tu = linspace(0,2*pi,n_uniform+1)';
    tu = tu(1:end-1) + pi/max(n_uniform,1);
    tu = mod(tu,2*pi);

    tc1 = local_cluster_angles(n_cluster_1);
    tc2 = local_cluster_angles(n_cluster_2);

    t_body = mod([tu; theta_1 + tc1; theta_2 + tc2],2*pi);
    t_body = sort(t_body);
    if numel(t_body) > nout
        t_body = t_body(1:nout);
    end
    rout_cells{body} = q_trip(body) + exp(1i*t_body);
end
end

function tc = local_cluster_angles(nc)
if nc <= 0
    tc = zeros(0,1);
    return
end
r_mob = 0.88;
npool = max(8*nc,128);
tseed = linspace(0,2*pi,npool+1)';
tseed = tseed(1:end-1);
tmob = sort(mobius_angle_map(tseed,r_mob));
dgap = min(abs(tmob),2*pi-abs(tmob));
[~,ord] = sort(dgap,'ascend');
tc = tmob(ord(1:nc));
end

function theta = mobius_angle_map(t,r)
z = exp(1i*t);
w = (z-r)./(1-r*z);
theta = mod(angle(w),2*pi);
end

function [abs_inf,abs_pw,rel_inf,rel_pw] = triplet_error_metrics(err_vec,ref_vec)
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
