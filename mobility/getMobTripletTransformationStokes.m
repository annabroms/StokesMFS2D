function [rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y, ...
    tau_stokes_nonpx,tau_stokes_nonpy,tau_stokes_e_nonpx,tau_stokes_e_nonpy, ...
    rimage_k,diag_info] = getMobTripletTransformationStokes(tau,geom,basis)
%GETMOBTRIPLETTRANSFORMATIONSTOKES Map mobility data to explicit 3-body sources.

rvec_out = geom.rvec_out;
q = geom.q;
triplets = geom.triplets;
rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
opt = geom.opt;
triplet_mode = lower(char(getOptField(opt,'triplet_accumulation_mode','edge_weighted')));

U = basis.U;
Y = basis.Y;
Lc = basis.Lc;
triplet_cache = basis.triplet_cache;
if isfield(geom,'triplet_weights') && ~isempty(geom.triplet_weights)
    triplet_weights = geom.triplet_weights;
else
    triplet_weights = ones(size(triplets,1),3);
end

P = numel(q);
N_c = opt.N_c;
N_f = opt.N_f;
N_large = numel(rvec_out)/P;
PM = numel(rvec_out);

n_coarse_tot = P*N_c;
tau_stokes_x_coarse = zeros(n_coarse_tot,1);
tau_stokes_y_coarse = zeros(n_coarse_tot,1);
tau_stokes_nonpx_coarse = zeros(n_coarse_tot,1);
tau_stokes_nonpy_coarse = zeros(n_coarse_tot,1);
rvec_in_coarse = zeros(n_coarse_tot,1);
coarse_ind = cell(P,1);

for i = 1:P
    tau_particle_x = tau((i-1)*N_large+1:i*N_large);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+i*N_large);

    step1 = U{1}*[tau_particle_x; tau_particle_y];
    lambda_nonproj_i = Y{1}*step1;
    lambda_proj_i = lambda_nonproj_i - Lc*lambda_nonproj_i;

    coarse_ind_i = (i-1)*N_c+1:i*N_c;
    coarse_ind{i} = coarse_ind_i;

    tau_stokes_nonpx_coarse(coarse_ind_i) = lambda_nonproj_i(1:N_c);
    tau_stokes_nonpy_coarse(coarse_ind_i) = lambda_nonproj_i(N_c+1:end);
    tau_stokes_x_coarse(coarse_ind_i) = lambda_proj_i(1:N_c);
    tau_stokes_y_coarse(coarse_ind_i) = lambda_proj_i(N_c+1:end);
    rvec_in_coarse(coarse_ind_i) = q(i)+rbase_in_c;
end

diag_info = init_diag_info(strcmp(triplet_mode,'pair_subtracted'));
acc = init_mobility_accumulator(P);

if strcmp(triplet_mode,'pair_subtracted')
    [~,~,tau_pair_x,tau_pair_y,tau_pair_nonpx,tau_pair_nonpy, ...
        tau_pair_e_nonpx,tau_pair_e_nonpy,rimage_pair] = ...
        getMobPairTransformationStokes(tau,geom,basis);
    acc = seed_pair_accumulator(acc,geom,tau_pair_x,tau_pair_y,tau_pair_nonpx, ...
        tau_pair_nonpy,tau_pair_e_nonpx,tau_pair_e_nonpy,rimage_pair);
    active_particles = sort(unique(geom.pairs(:)));
else
    active_particles = sort(unique(triplets(:)));
end

for row = 1:size(triplets,1)
    if strcmp(triplet_mode,'pair_subtracted')
        [triplet,idx,beta_trip,diag_info] = build_pair_subtracted_triplet_beta( ...
            row,geom,basis,coarse_ind,tau_stokes_x_coarse,tau_stokes_y_coarse,diag_info);
        body_weights = ones(1,3);
    else
        triplet = getStokesTripletInstance(triplet_cache,row);
        im_counts = cellfun(@numel,triplet.rimage_cells(:)).';
        idx = get_triplet_indices(N_f,im_counts);
        rhs_trip = pack_triplet_coarse_rhs(triplet.meta.ids(:),coarse_ind, ...
            tau_stokes_x_coarse,tau_stokes_y_coarse);
        rhs_trip = rotatePairOrderedStokesData(rhs_trip,N_c, ...
            triplet.meta.phase_c,conj(triplet.meta.rot));
        beta_ref = triplet.group.Ypf*(triplet.group.Upf*rhs_trip);
        beta_trip = rotateStokesTripletSourceVector(beta_ref,N_f,im_counts, ...
            triplet.meta.phase_f_inv,triplet.meta.rot);
        body_weights = triplet_weights(row,:);
    end

    for body = 1:3
        pid = triplet.meta.ids(body);
        weight = body_weights(body);
        rsrc_body = [q(pid)+rbase_in_f; triplet.rimage_cells{body}];
        body_nonproj = extract_triplet_body_sources(beta_trip,idx,body);
        body_nonproj = scale_body_sources(body_nonproj,weight);
        tau_nonproj = pack_body_sources(body_nonproj);
        tau_proj = projectOutRigid2D(tau_nonproj,rsrc_body,q(pid));
        body_proj = split_projected_body(tau_proj,N_f,numel(triplet.rimage_cells{body}));

        acc.fine_x{pid} = add_or_init(acc.fine_x{pid},body_proj.fine_x);
        acc.fine_y{pid} = add_or_init(acc.fine_y{pid},body_proj.fine_y);
        acc.fine_nonpx{pid} = add_or_init(acc.fine_nonpx{pid},body_nonproj.fine_x);
        acc.fine_nonpy{pid} = add_or_init(acc.fine_nonpy{pid},body_nonproj.fine_y);

        acc.extra_x_chunks{pid}{end+1,1} = body_proj.extra_x;
        acc.extra_y_chunks{pid}{end+1,1} = body_proj.extra_y;
        acc.extra_nonpx_chunks{pid}{end+1,1} = body_nonproj.extra_x;
        acc.extra_nonpy_chunks{pid}{end+1,1} = body_nonproj.extra_y;
        acc.rimage_chunks{pid}{end+1,1} = triplet.rimage_cells{body};
    end
end

[rvec_in,tau_stokes_x,tau_stokes_y,tau_stokes_nonpx,tau_stokes_nonpy, ...
    tau_stokes_e_nonpx,tau_stokes_e_nonpy,rimage_k] = finalize_mobility_accumulator( ...
    acc,active_particles,q,rbase_in_c,rbase_in_f, ...
    tau_stokes_x_coarse,tau_stokes_y_coarse, ...
    tau_stokes_nonpx_coarse,tau_stokes_nonpy_coarse,rvec_in_coarse,N_f);

end

function acc = init_mobility_accumulator(P)
acc = struct();
acc.fine_x = cell(P,1);
acc.fine_y = cell(P,1);
acc.fine_nonpx = cell(P,1);
acc.fine_nonpy = cell(P,1);
acc.extra_x_chunks = repmat({cell(0,1)},P,1);
acc.extra_y_chunks = repmat({cell(0,1)},P,1);
acc.extra_nonpx_chunks = repmat({cell(0,1)},P,1);
acc.extra_nonpy_chunks = repmat({cell(0,1)},P,1);
acc.rimage_chunks = repmat({cell(0,1)},P,1);
end

function acc = seed_pair_accumulator(acc,geom,tau_x,tau_y,tau_nonpx,tau_nonpy, ...
    tau_e_nonpx,tau_e_nonpy,rimage_k)
P = numel(geom.q);
N_c = geom.opt.N_c;
N_f = geom.opt.N_f;
active_pairs = sort(unique(geom.pairs(:)));
offset_xy = P*N_c;
offset_nonp = P*N_c;

for it = 1:numel(active_pairs)
    pid = active_pairs(it);

    fine_xy_idx = offset_xy + (1:N_f);
    acc.fine_x{pid} = tau_x(fine_xy_idx);
    acc.fine_y{pid} = tau_y(fine_xy_idx);
    offset_xy = offset_xy + N_f;

    fine_nonp_idx = offset_nonp + (1:N_f);
    acc.fine_nonpx{pid} = tau_nonpx(fine_nonp_idx);
    acc.fine_nonpy{pid} = tau_nonpy(fine_nonp_idx);
    offset_nonp = offset_nonp + N_f;

    if ~isempty(rimage_k{pid})
        n_extra = numel(rimage_k{pid});
        extra_xy_idx = offset_xy + (1:n_extra);
        acc.extra_x_chunks{pid}{1,1} = tau_x(extra_xy_idx);
        acc.extra_y_chunks{pid}{1,1} = tau_y(extra_xy_idx);
        acc.extra_nonpx_chunks{pid}{1,1} = tau_e_nonpx{pid};
        acc.extra_nonpy_chunks{pid}{1,1} = tau_e_nonpy{pid};
        acc.rimage_chunks{pid}{1,1} = rimage_k{pid};
        offset_xy = offset_xy + n_extra;
    end
end
end

function [triplet,idx,beta_3only,diag_info] = build_pair_subtracted_triplet_beta( ...
    row,geom,basis,coarse_ind,tau_x_coarse,tau_y_coarse,diag_info)
triplet = getStokesTripletInstance(basis.triplet_cache,row);
ids = triplet.meta.ids(:);
N_c = geom.opt.N_c;
N_f = geom.opt.N_f;
im_counts = cellfun(@numel,triplet.rimage_cells(:)).';
idx = get_triplet_indices(N_f,im_counts);

rhs_trip = pack_triplet_coarse_rhs(ids,coarse_ind,tau_x_coarse,tau_y_coarse);
rhs_trip = rotatePairOrderedStokesData(rhs_trip,N_c,triplet.meta.phase_c,conj(triplet.meta.rot));
if isfield(triplet.group,'Upf_full') && ~isempty(triplet.group.Upf_full)
    beta_ref = triplet.group.Ypf_full*(triplet.group.Upf_full*rhs_trip);
else
    beta_ref = triplet.group.Ypf*(triplet.group.Upf*rhs_trip);
end
beta_full = rotateStokesTripletSourceVector(beta_ref,N_f,im_counts, ...
    triplet.meta.phase_f_inv,triplet.meta.rot);

beta_pair_sum = zeros(size(beta_full));
edge_list = [1 2; 1 3; 2 3];
for edge = 1:size(edge_list,1)
    local_a = edge_list(edge,1);
    local_b = edge_list(edge,2);
    pair_row = geom.pair_lookup(ids(local_a),ids(local_b));
    if pair_row == 0
        error('Missing pair row for triplet edge (%d,%d).',ids(local_a),ids(local_b));
    end

    beta_pair = get_pair_beta(pair_row,geom,basis,coarse_ind,tau_x_coarse,tau_y_coarse);
    beta_pair_sum = beta_pair_sum + embed_pair_into_triplet(beta_pair,pair_row, ...
        triplet,idx,geom);
end

beta_3only = beta_full - beta_pair_sum;
diag_info = update_diag_info(diag_info,row,ids,beta_full,beta_pair_sum,beta_3only);
end

function beta_pair = get_pair_beta(pair_row,geom,basis,coarse_ind,tau_x_coarse,tau_y_coarse)
pair_ids = geom.pairs(pair_row,:);
N_c = geom.opt.N_c;
N_f = geom.opt.N_f;

lambda_i = [tau_x_coarse(coarse_ind{pair_ids(1)}); tau_y_coarse(coarse_ind{pair_ids(1)})];
lambda_j = [tau_x_coarse(coarse_ind{pair_ids(2)}); tau_y_coarse(coarse_ind{pair_ids(2)})];

if isfield(basis,'pair_cache') && isfield(basis.pair_cache,'enabled') && basis.pair_cache.enabled
    pair = getStokesPairInstance(basis.pair_cache,pair_row);
    rhs_pair = [lambda_i(1:N_c) lambda_j(1:N_c) lambda_i(N_c+1:end) lambda_j(N_c+1:end)];
    rhs_pair = rotatePairOrderedStokesData(rhs_pair,N_c,pair.meta.phase_c,conj(pair.meta.rot));
    beta_ref = pair.group.Ypf*(pair.group.Upf*rhs_pair(:));
    beta_pair = rotateStokesPairSourceVector(beta_ref,N_f,numel(pair.rimage_i), ...
        numel(pair.rimage_j),pair.meta.phase_f_inv,pair.meta.rot);
    return
end

if isempty(basis.Upf) || isempty(basis.Ypf)
    error(['Pair-subtracted triplets require either pair_cache or explicit ', ...
        'pair Upf/Ypf data in basis.']);
end

rhs_pair = [lambda_i(1:N_c); lambda_j(1:N_c); lambda_i(N_c+1:end); lambda_j(N_c+1:end)];
beta_pair = basis.Ypf{pair_ids(1),pair_ids(2)} * (basis.Upf{pair_ids(1),pair_ids(2)} * rhs_pair);
end

function beta_embed = embed_pair_into_triplet(beta_pair,pair_row,triplet,idx,geom)
pair_ids = geom.pairs(pair_row,:);
im_i = numel(geom.rimage_vec{pair_ids(1),pair_ids(2)});
im_j = numel(geom.rimage_vec{pair_ids(2),pair_ids(1)});
pair_idx = get_pair_indices(geom.opt.N_f,im_i,im_j);
pair_body_1 = extract_pair_body_sources(beta_pair,pair_idx,1);
pair_body_2 = extract_pair_body_sources(beta_pair,pair_idx,2);

beta_embed = zeros(idx.ntot,1);
beta_embed = assign_pair_body(beta_embed,triplet,idx,geom,pair_ids(1),pair_ids(2),pair_body_1);
beta_embed = assign_pair_body(beta_embed,triplet,idx,geom,pair_ids(2),pair_ids(1),pair_body_2);
end

function beta_embed = assign_pair_body(beta_embed,triplet,idx,geom,body_id,neigh_id,body_src)
local_body = find(triplet.meta.ids == body_id,1,'first');
if isempty(local_body)
    error('Triplet body %d not found during pair embedding.',body_id);
end

[extra_idx_x,extra_idx_y] = get_triplet_edge_chunk_indices(triplet,idx,geom,local_body,neigh_id);
[fine_idx_x,fine_idx_y] = get_triplet_fine_indices(idx,local_body);

beta_embed(fine_idx_x) = beta_embed(fine_idx_x) + body_src.fine_x;
beta_embed(extra_idx_x) = beta_embed(extra_idx_x) + body_src.extra_x;
beta_embed(fine_idx_y) = beta_embed(fine_idx_y) + body_src.fine_y;
beta_embed(extra_idx_y) = beta_embed(extra_idx_y) + body_src.extra_y;
end

function [extra_idx_x,extra_idx_y] = get_triplet_edge_chunk_indices(triplet,idx,geom,local_body,neigh_id)
ids = triplet.meta.ids(:).';
neigh_local = setdiff(1:3,local_body,'stable');
neigh_ids = ids(neigh_local);
counts = [numel(geom.rimage_vec{ids(local_body),neigh_ids(1)}), ...
    numel(geom.rimage_vec{ids(local_body),neigh_ids(2)})];

if neigh_id == neigh_ids(1)
    offset = 0;
    count = counts(1);
elseif neigh_id == neigh_ids(2)
    offset = counts(1);
    count = counts(2);
else
    error('Neighbour %d not found in triplet chunk lookup.',neigh_id);
end

[body_extra_x,body_extra_y] = get_triplet_body_extra_indices(idx,local_body);
if count == 0
    extra_idx_x = zeros(0,1);
    extra_idx_y = zeros(0,1);
else
    extra_idx_x = body_extra_x(offset+1:offset+count);
    extra_idx_y = body_extra_y(offset+1:offset+count);
end
end

function [fine_idx_x,fine_idx_y] = get_triplet_fine_indices(idx,body)
switch body
    case 1
        fine_idx_x = idx.f1x;
        fine_idx_y = idx.f1y;
    case 2
        fine_idx_x = idx.f2x;
        fine_idx_y = idx.f2y;
    case 3
        fine_idx_x = idx.f3x;
        fine_idx_y = idx.f3y;
    otherwise
        error('Invalid triplet body index.');
end
end

function [extra_idx_x,extra_idx_y] = get_triplet_body_extra_indices(idx,body)
switch body
    case 1
        extra_idx_x = idx.e1x;
        extra_idx_y = idx.e1y;
    case 2
        extra_idx_x = idx.e2x;
        extra_idx_y = idx.e2y;
    case 3
        extra_idx_x = idx.e3x;
        extra_idx_y = idx.e3y;
    otherwise
        error('Invalid triplet body index.');
end
end

function body_src = extract_pair_body_sources(beta_pair,idx,body)
switch body
    case 1
        body_src.fine_x = beta_pair(idx.f1x);
        body_src.extra_x = beta_pair(idx.e1x);
        body_src.fine_y = beta_pair(idx.f1y);
        body_src.extra_y = beta_pair(idx.e1y);
    case 2
        body_src.fine_x = beta_pair(idx.f2x);
        body_src.extra_x = beta_pair(idx.e2x);
        body_src.fine_y = beta_pair(idx.f2y);
        body_src.extra_y = beta_pair(idx.e2y);
    otherwise
        error('Invalid pair body index.');
end
end

function body_src = extract_triplet_body_sources(beta,idx,body)
switch body
    case 1
        body_src.fine_x = beta(idx.f1x);
        body_src.extra_x = beta(idx.e1x);
        body_src.fine_y = beta(idx.f1y);
        body_src.extra_y = beta(idx.e1y);
    case 2
        body_src.fine_x = beta(idx.f2x);
        body_src.extra_x = beta(idx.e2x);
        body_src.fine_y = beta(idx.f2y);
        body_src.extra_y = beta(idx.e2y);
    case 3
        body_src.fine_x = beta(idx.f3x);
        body_src.extra_x = beta(idx.e3x);
        body_src.fine_y = beta(idx.f3y);
        body_src.extra_y = beta(idx.e3y);
    otherwise
        error('Invalid triplet body index.');
end
end

function body_src = scale_body_sources(body_src,weight)
body_src.fine_x = weight*body_src.fine_x;
body_src.extra_x = weight*body_src.extra_x;
body_src.fine_y = weight*body_src.fine_y;
body_src.extra_y = weight*body_src.extra_y;
end

function tau_body = pack_body_sources(body_src)
tau_body = [body_src.fine_x; body_src.extra_x; body_src.fine_y; body_src.extra_y];
end

function body_proj = split_projected_body(tau_proj,Nf,nim)
body_proj = struct();
body_proj.fine_x = tau_proj(1:Nf);
body_proj.extra_x = tau_proj(Nf+1:Nf+nim);
body_proj.fine_y = tau_proj(Nf+nim+1:2*Nf+nim);
body_proj.extra_y = tau_proj(2*Nf+nim+1:2*Nf+2*nim);
end

function rhs_trip = pack_triplet_coarse_rhs(ids,coarse_ind,tau_x,tau_y)
rhs_trip = [];
for body = 1:3
    rhs_trip = [rhs_trip; tau_x(coarse_ind{ids(body)})];
end
for body = 1:3
    rhs_trip = [rhs_trip; tau_y(coarse_ind{ids(body)})];
end
end

function idx = get_triplet_indices(Nf,im_counts)
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

function idx = get_pair_indices(Nf,im_i,im_j)
idx = struct();
idx.f1x = 1:Nf;
idx.e1x = Nf+1:Nf+im_i;
idx.f2x = Nf+im_i+1:2*Nf+im_i;
idx.e2x = 2*Nf+im_i+1:2*Nf+im_i+im_j;
idx.f1y = 2*Nf+im_i+im_j+1:3*Nf+im_i+im_j;
idx.e1y = 3*Nf+im_i+im_j+1:3*Nf+2*im_i+im_j;
idx.f2y = 3*Nf+2*im_i+im_j+1:4*Nf+2*im_i+im_j;
idx.e2y = 4*Nf+2*im_i+im_j+1:4*Nf+2*im_i+2*im_j;
end

function [rvec_in,tau_x,tau_y,tau_nonpx,tau_nonpy,tau_e_nonpx,tau_e_nonpy,rimage_k] = ...
    finalize_mobility_accumulator(acc,active_particles,q,rbase_in_c,rbase_in_f, ...
    tau_x_coarse,tau_y_coarse,tau_nonpx_coarse,tau_nonpy_coarse,rvec_in_coarse,N_f)
[rimage_k,tau_e_x,tau_e_y,tau_e_nonpx,tau_e_nonpy] = collect_extra_sources(acc);

x_chunks = cell(1+2*numel(active_particles),1);
y_chunks = cell(1+2*numel(active_particles),1);
nonpx_chunks = cell(1+numel(active_particles),1);
nonpy_chunks = cell(1+numel(active_particles),1);
rvec_chunks = cell(1+2*numel(active_particles),1);

x_chunks{1} = tau_x_coarse;
y_chunks{1} = tau_y_coarse;
nonpx_chunks{1} = tau_nonpx_coarse;
nonpy_chunks{1} = tau_nonpy_coarse;
rvec_chunks{1} = rvec_in_coarse;

for it = 1:numel(active_particles)
    pid = active_particles(it);
    x_chunks{2*it} = ensure_fine_column(acc.fine_x{pid},N_f);
    y_chunks{2*it} = ensure_fine_column(acc.fine_y{pid},N_f);
    nonpx_chunks{it+1} = ensure_fine_column(acc.fine_nonpx{pid},N_f);
    nonpy_chunks{it+1} = ensure_fine_column(acc.fine_nonpy{pid},N_f);
    rvec_chunks{2*it} = q(pid)+rbase_in_f;
    x_chunks{2*it+1} = tau_e_x{pid};
    y_chunks{2*it+1} = tau_e_y{pid};
    rvec_chunks{2*it+1} = rimage_k{pid};
end

tau_x = vertcat(x_chunks{:});
tau_y = vertcat(y_chunks{:});
tau_nonpx = vertcat(nonpx_chunks{:});
tau_nonpy = vertcat(nonpy_chunks{:});
rvec_in = vertcat(rvec_chunks{:});
end

function [rimage_k,ex,ey,enx,eny] = collect_extra_sources(acc)
P = numel(acc.fine_x);
rimage_k = cell(P,1);
ex = cell(P,1);
ey = cell(P,1);
enx = cell(P,1);
eny = cell(P,1);
for k = 1:P
    if isempty(acc.rimage_chunks{k})
        rimage_k{k} = zeros(0,1);
        ex{k} = zeros(0,1);
        ey{k} = zeros(0,1);
        enx{k} = zeros(0,1);
        eny{k} = zeros(0,1);
    else
        rimage_k{k} = vertcat(acc.rimage_chunks{k}{:});
        ex{k} = vertcat(acc.extra_x_chunks{k}{:});
        ey{k} = vertcat(acc.extra_y_chunks{k}{:});
        enx{k} = vertcat(acc.extra_nonpx_chunks{k}{:});
        eny{k} = vertcat(acc.extra_nonpy_chunks{k}{:});
    end
end
end

function vec = add_or_init(vec,contrib)
if isempty(vec)
    vec = contrib;
else
    vec = vec + contrib;
end
end

function vec = ensure_fine_column(vec,Nf)
if isempty(vec)
    vec = zeros(Nf,1);
else
    vec = vec(:);
end
end

function diag_info = init_diag_info(use_pair_subtraction)
diag_info = struct();
diag_info.used_pair_subtraction = use_pair_subtraction;
diag_info.rep_row = 0;
diag_info.rep_ids = zeros(1,3);
diag_info.rep_norm_full = 0;
diag_info.rep_norm_pair_sum = 0;
diag_info.rep_norm_remainder = 0;
diag_info.max_norm_full = 0;
diag_info.max_norm_pair_sum = 0;
diag_info.max_norm_remainder = 0;
end

function diag_info = update_diag_info(diag_info,row,ids,beta_full,beta_pair_sum,beta_3only)
nf = norm(beta_full);
np = norm(beta_pair_sum);
nr = norm(beta_3only);

if diag_info.rep_row == 0
    diag_info.rep_row = row;
    diag_info.rep_ids = ids(:).';
    diag_info.rep_norm_full = nf;
    diag_info.rep_norm_pair_sum = np;
    diag_info.rep_norm_remainder = nr;
end

diag_info.max_norm_full = max(diag_info.max_norm_full,nf);
diag_info.max_norm_pair_sum = max(diag_info.max_norm_pair_sum,np);
diag_info.max_norm_remainder = max(diag_info.max_norm_remainder,nr);
end
