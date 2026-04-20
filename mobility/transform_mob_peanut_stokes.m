function [lam_c_x, lam_c_nonpx,lam_self_x, lam_f_x, ...
    lam_c_y,lam_c_nonpy,lam_self_y,lam_f_y,u_corr,rimage_k] = ...
    transform_mob_peanut_stokes(tau,geom,basis)
%TRANSFORM_MOB_PEANUT_STOKES Map boundary data to coarse interior Stokes sources.
%
% Syntax:
%   [lam_c_x, lam_c_nonpx, lam_self_x, lam_f_x, ...
%    lam_c_y, lam_c_nonpy, lam_self_y, lam_f_y, u_corr, rimage_k] = ...
%       transform_mob_peanut_stokes(tau,geom,basis)
%
% Outputs:
%   lam_c_x, lam_c_y
%       Projected coarse source strengths with compressed pair corrections.
%   lam_c_nonpx, lam_c_nonpy
%       Unprojected coarse source strengths, needed for rigid-body terms 
%   lam_self_x, lam_self_y
%       Coarse projected sources before pair corrections are applied
%       (corresponding to one-body basis only). Needed to enforce identity
%       diagonal blocks in system matrix.
%   lam_f_x, lam_f_y
%       Fine source strengths used in the pair correction and, when
%       `opt.cmap = 0`, also in the BK' correction path.
%   u_corr
%       Pair correction to the velocity field evaluated at all collocation nodes. 
%   rimage_k
%       Cell array of per-body enhancing source locations.
%
%   Anna Broms, Mar 2026
%

%% Prepare discretisation and preallocate
rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
rimage_vec = geom.rimage_vec;
opt = geom.opt;
rvec_out = geom.rvec_out;
rcheck_out = geom.rcheck;
q = geom.q;
pairs = geom.pairs;

U = basis.U;
Y = basis.Y;
Lc = basis.Lc;
Upf = basis.Upf;
Ypf = basis.Ypf;
DC_all = basis.DC_all;
YC_all = basis.YC_all;
Cmap = basis.Cmap;
Lc_pair = basis.Lc_pair;
if isfield(basis,'pair_cache')
    pair_cache = basis.pair_cache;
else
    pair_cache = struct('enabled',false);
end
use_pair_cache = isfield(pair_cache,'enabled') && pair_cache.enabled;

P = length(q);
n_pairs = size(pairs,1);
N_c = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P;
PM = length(rvec_out);
N_check = length(rcheck_out)/P;
PM2 = length(rcheck_out);
use_cmap = logical(getOptField(opt,'cmap',false));
use_dense = logical(getOptField(opt,'use_dense',false)) && isequal(rcheck_out,rvec_out);
use_matrix_free_Lc_pair = logical(getOptField(opt,'use_matrix_free_Lc_pair',true));
get_bndry_field = logical(getOptField(opt,'get_bndry_field',true));
need_explicit_pair_sources = ~use_cmap || get_bndry_field || ~opt.self_correct;
has_pair_meta = isfield(pair_cache,'meta') && numel(pair_cache.meta) >= n_pairs;
if use_matrix_free_Lc_pair
    coarse_rigid_projector = buildRigidProjectionCache2D(rbase_in_c);
else
    coarse_rigid_projector = [];
end

lam_c_x = zeros(N_c*P,1);
lam_c_y = zeros(N_c*P,1);
lam_c_nonpx = zeros(N_c*P,1);
lam_c_nonpy = zeros(N_c*P,1);

if need_explicit_pair_sources
    lam_beta_f_x = zeros(N_f*P,1);
    lam_beta_f_y = zeros(N_f*P,1);
    lam_beta_e_x_chunks = repmat({cell(0,1)},P,1);
    lam_beta_e_y_chunks = repmat({cell(0,1)},P,1);
else
    lam_beta_f_x = [];
    lam_beta_f_y = [];
    lam_beta_e_x_chunks = {};
    lam_beta_e_y_chunks = {};
end

rimage_k = cell(P,1);
u_corr = zeros(2*N_check*P,1);

%% Phase 1: one-body coarse map.
for i = 1:P
    coarse_ind = (i-1)*N_c+1:i*N_c;
    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    step1 = U{1}*[tau_particle_x;tau_particle_y];
    tau_mapped = Y{1}*step1;

    lam_c_nonpx(coarse_ind) = tau_mapped(1:N_c);
    lam_c_nonpy(coarse_ind) = tau_mapped(N_c+1:end);

    tau_i_x = tau_mapped(1:N_c);
    tau_i_y = tau_mapped(N_c+1:end);
    tau_mapped = [tau_i_x; tau_i_y]-Lc*[tau_i_x; tau_i_y];

    lam_c_x(coarse_ind) = tau_mapped(1:N_c);
    lam_c_y(coarse_ind) = tau_mapped(N_c+1:end);
end

lam_self_x = lam_c_x;
lam_self_y = lam_c_y;


%% Phase 2: loop over all pairs
for row = 1:n_pairs
    i = pairs(row,1);
    p2 = pairs(row,2);
    if has_pair_meta
        meta = pair_cache.meta(row);
    else
        meta = struct();
    end

    coarse_i = (i-1)*N_c+1:i*N_c;
    coarse_p2 = (p2-1)*N_c+1:p2*N_c;

    % take out projected one-body sources
    lambda_i = [lam_self_x(coarse_i); lam_self_y(coarse_i)];
    lambda_p2 = [lam_self_x(coarse_p2); lam_self_y(coarse_p2)];

    rhs_mat = [lambda_i(1:N_c) lambda_p2(1:N_c) ...
               lambda_i(N_c+1:end) lambda_p2(N_c+1:end)];

    if use_pair_cache % reuse information from multiple pairs
        group = pair_cache.groups(meta.group_id);
        rhs = rotatePairOrderedStokesData(rhs_mat,N_c,meta.phase_c,conj(meta.rot));
        rhs = rhs(:);
        rimage_i = mapCanonicalPointsToActual(group.rimage_canon{1},meta);
        rimage_p2 = mapCanonicalPointsToActual(group.rimage_canon{2},meta);
    else
        rhs = rhs_mat(:);
        rimage_i = rimage_vec{i,p2};
        rimage_p2 = rimage_vec{p2,i};
    end

    rimage_k{i} = [rimage_k{i}; rimage_i];
    rimage_k{p2} = [rimage_k{p2}; rimage_p2];

    tau_mapped_tot = [];
    if need_explicit_pair_sources
        if use_pair_cache
            pair_mapped = group.Upf*rhs;
            tau_mapped_loc = group.Ypf*pair_mapped;
            tau_mapped_tot = rotateStokesPairSourceVector(tau_mapped_loc,N_f,length(rimage_i),length(rimage_p2), ...
                meta.phase_f_inv,meta.rot);
        else
            pair_mapped = Upf{i,p2}*rhs;
            tau_mapped_tot = Ypf{i,p2}*pair_mapped;
            tau_mapped_loc = tau_mapped_tot;
        end
    else
        tau_mapped_loc = [];
    end

    if use_cmap
        if use_pair_cache
            tau_peanut_ntot = group.Cmap*rhs;
        else
            tau_peanut_ntot = Cmap{i,p2}*rhs;
        end
    else
        if use_pair_cache
            tau_peanut_ntot = group.YC*(group.DC*tau_mapped_loc);
        else
            tau_peanut_ntot = YC_all{i,p2}*(DC_all{i,p2}*tau_mapped_tot);
        end
    end

    if use_matrix_free_Lc_pair
        if use_pair_cache
            tau_peanut_loc = projectOutRigidPairCached2D( ...
                tau_peanut_ntot,coarse_rigid_projector);
            tau_peanut_tot = rotatePairOrderedStokesData( ...
                tau_peanut_loc,N_c,meta.phase_c_inv,meta.rot);
            tau_peanut_tot = tau_peanut_tot(:);
        else
            tau_peanut_tot = projectOutRigidPairCached2D( ...
                tau_peanut_ntot,coarse_rigid_projector);
        end
    else
        if use_pair_cache
            tau_peanut_loc = group.Lc_pair*tau_peanut_ntot;
            tau_peanut_tot = rotatePairOrderedStokesData( ...
                tau_peanut_loc,N_c,meta.phase_c_inv,meta.rot);
            tau_peanut_tot = tau_peanut_tot(:);
        else
            tau_peanut_tot = Lc_pair*tau_peanut_ntot;
        end
    end

    lam_c_x(coarse_i) = lam_c_x(coarse_i) + tau_peanut_tot(1:N_c);
    lam_c_x(coarse_p2) = lam_c_x(coarse_p2) + tau_peanut_tot(N_c+1:2*N_c);
    lam_c_y(coarse_i) = lam_c_y(coarse_i) + tau_peanut_tot(2*N_c+1:3*N_c);
    lam_c_y(coarse_p2) = lam_c_y(coarse_p2) + tau_peanut_tot(3*N_c+1:4*N_c);

    if need_explicit_pair_sources
        im_nr = length(rimage_i);
        f_ind1_x = 1:N_f;
        e_ind1_x = N_f+1:N_f+im_nr;
        f_ind2_x = N_f+im_nr+1:2*N_f+im_nr;
        e_ind2_x = 2*N_f+im_nr+1:2*N_f+2*im_nr;
        f_ind1_y = 2*N_f+2*im_nr+1:3*N_f+2*im_nr;
        e_ind1_y = 3*N_f+2*im_nr+1:3*N_f+3*im_nr;
        f_ind2_y = 3*N_f+3*im_nr+1:4*N_f+3*im_nr;
        e_ind2_y = 4*N_f+3*im_nr+1:4*N_f+4*im_nr;

        fine_i = (i-1)*N_f+1:i*N_f;
        fine_p2 = (p2-1)*N_f+1:p2*N_f;
        lam_beta_f_x(fine_i) = lam_beta_f_x(fine_i) + tau_mapped_tot(f_ind1_x);
        lam_beta_f_x(fine_p2) = lam_beta_f_x(fine_p2) + tau_mapped_tot(f_ind2_x);
        lam_beta_f_y(fine_i) = lam_beta_f_y(fine_i) + tau_mapped_tot(f_ind1_y);
        lam_beta_f_y(fine_p2) = lam_beta_f_y(fine_p2) + tau_mapped_tot(f_ind2_y);

        lam_beta_e_x_chunks{i}{end+1,1} = tau_mapped_tot(e_ind1_x);
        lam_beta_e_x_chunks{p2}{end+1,1} = tau_mapped_tot(e_ind2_x);
        lam_beta_e_y_chunks{i}{end+1,1} = tau_mapped_tot(e_ind1_y);
        lam_beta_e_y_chunks{p2}{end+1,1} = tau_mapped_tot(e_ind2_y);
    end

    rcheck_i = rcheck_out((i-1)*N_check+1:i*N_check,:);
    rcheck_p2 = rcheck_out((p2-1)*N_check+1:p2*N_check,:);
    rout_pair = [rcheck_i; rcheck_p2];
    rin_pair_c = [rbase_in_c+q(i); rbase_in_c+q(p2)];

    if  ~opt.self_correct
        % When the fine pair sources are available, evaluate the pair
        % correction directly from the projected fine densities to mirror
        % transform_mob_peanut.
        fine_rin_pair = [rbase_in_f+q(i); rimage_i; rbase_in_f+q(p2); rimage_p2];

        nsrc_i = numel(rbase_in_f) + numel(rimage_i);
        nsrc_p2 = numel(rbase_in_f) + numel(rimage_p2);

        tau_fine_i = [tau_mapped_tot(f_ind1_x); tau_mapped_tot(e_ind1_x); ...
                      tau_mapped_tot(f_ind1_y); tau_mapped_tot(e_ind1_y)];
        tau_fine_p2 = [tau_mapped_tot(f_ind2_x); tau_mapped_tot(e_ind2_x); ...
                       tau_mapped_tot(f_ind2_y); tau_mapped_tot(e_ind2_y)];

        tau_mapped_proj_i = projectOutRigid2D(tau_fine_i, ...
            [rbase_in_f+q(i); rimage_i], q(i));
        tau_mapped_proj_p2 = projectOutRigid2D(tau_fine_p2, ...
            [rbase_in_f+q(p2); rimage_p2], q(p2));

        fine_fx = [tau_mapped_proj_i(1:nsrc_i); tau_mapped_proj_p2(1:nsrc_p2)];
        fine_fy = [tau_mapped_proj_i(nsrc_i+1:end); tau_mapped_proj_p2(nsrc_p2+1:end)];
        if use_dense && isfield(meta,'Upair_colloc_actual') && ~isempty(meta.Upair_colloc_actual)
            u_pair = meta.Upair_colloc_actual * [fine_fx; fine_fy];
        else
            [u1,v1] = stokSLPdirect(real(fine_rin_pair),imag(fine_rin_pair), ...
                real(rout_pair),imag(rout_pair),fine_fx,fine_fy,numel(fine_rin_pair));
            u_pair = [u1; v1];
        end
    else

        rhs_self = [lam_self_x(coarse_i); lam_self_x(coarse_p2); ...
                    lam_self_y(coarse_i); lam_self_y(coarse_p2)];
        if use_dense && isfield(meta,'Ucross_colloc_actual') && ~isempty(meta.Ucross_colloc_actual)
            u_pair = meta.Ucross_colloc_actual * rhs_self;
        else
            [ui,vi] = stokSLPdirect(real(rbase_in_c+q(p2)),imag(rbase_in_c+q(p2)), ...
                real(rcheck_i),imag(rcheck_i),lam_self_x(coarse_p2),lam_self_y(coarse_p2),N_c);
            [up2,vp2] = stokSLPdirect(real(rbase_in_c+q(i)),imag(rbase_in_c+q(i)), ...
                real(rcheck_p2),imag(rcheck_p2),lam_self_x(coarse_i),lam_self_y(coarse_i),N_c);
            u_pair = -[ui; up2; vi; vp2];
        end
        % need to remove lr in the matvec to create identity in system matrix

    end

    if use_dense && isfield(meta,'Ecolloc_actual') && ~isempty(meta.Ecolloc_actual)
        tau_peanut_src = [tau_peanut_tot(1:2*N_c); tau_peanut_tot(2*N_c+1:4*N_c)];
        u_peanut_corr = meta.Ecolloc_actual * tau_peanut_src;
    else
        [u1,v1] = stokSLPdirect(real(rin_pair_c),imag(rin_pair_c), ...
            real(rout_pair),imag(rout_pair),tau_peanut_tot(1:2*N_c), ...
            tau_peanut_tot(2*N_c+1:4*N_c),2*N_c);
        u_peanut_corr = [u1; v1];
    end

    pair_ind = [(i-1)*N_check+1:i*N_check (p2-1)*N_check+1:p2*N_check ...
        (i-1)*N_check+PM2+1:i*N_check+PM2 ...
        (p2-1)*N_check+PM2+1:p2*N_check+PM2]';
    u_corr(pair_ind) = u_corr(pair_ind)+u_pair-u_peanut_corr;
end


%% Store fine sources for post-processing? 
if need_explicit_pair_sources
    lam_f_x = cell(P,1);
    lam_f_y = cell(P,1);
    for k = 1:P
        fine_ind = (k-1)*N_f+1:k*N_f;
        if isempty(lam_beta_e_x_chunks{k})
            lam_f_x{k} = lam_beta_f_x(fine_ind);
            lam_f_y{k} = lam_beta_f_y(fine_ind);
        else
            lam_f_x{k} = [lam_beta_f_x(fine_ind); vertcat(lam_beta_e_x_chunks{k}{:})];
            lam_f_y{k} = [lam_beta_f_y(fine_ind); vertcat(lam_beta_e_y_chunks{k}{:})];
        end
    end
else
    lam_f_x = [];
    lam_f_y = [];
end

end

function geom = buildRigidProjectionCache2D(rbase)
rbase = rbase(:);
dx = real(rbase);
dy = imag(rbase);
n = numel(rbase);

g13 = -sum(dy);
g23 = sum(dx);
g33 = sum(dx.^2 + dy.^2);
G = [n, 0, g13; ...
     0, n, g23; ...
     g13, g23, g33];

geom.n = n;
geom.dx = dx;
geom.dy = dy;
geom.Ginv = G\eye(3);
end

function pair_proj = projectOutRigidPairCached2D(pair_vec,geom)
n = geom.n;
dx = geom.dx;
dy = geom.dy;

x_i = pair_vec(1:n);
x_j = pair_vec(n+1:2*n);
y_i = pair_vec(2*n+1:3*n);
y_j = pair_vec(3*n+1:4*n);

alpha = geom.Ginv * [sum(x_i), sum(x_j); ...
                     sum(y_i), sum(y_j); ...
                     -dy.'*x_i + dx.'*y_i, -dy.'*x_j + dx.'*y_j];

pair_proj = pair_vec;
pair_proj(1:n) = x_i - alpha(1,1) + alpha(3,1)*dy;
pair_proj(n+1:2*n) = x_j - alpha(1,2) + alpha(3,2)*dy;
pair_proj(2*n+1:3*n) = y_i - alpha(2,1) - alpha(3,1)*dx;
pair_proj(3*n+1:4*n) = y_j - alpha(2,2) - alpha(3,2)*dx;
end

function z_actual = mapCanonicalPointsToActual(z_canon,meta)
z_canon = z_canon(:);
if isempty(z_canon)
    z_actual = zeros(0,1);
    return
end

z_actual = meta.mid + meta.rot*z_canon;
end
