function [lam_c_x, lam_c_nonpx, lam_self_x, lam_f_x, ...
    lam_c_y, lam_c_nonpy, lam_self_y, lam_f_y, u_corr, rimage_k] = ...
    transform_mob_peanut_big_sparse_stokes(tau,geom,basis)
%TRANSFORM_MOB_PEANUT_BIG_SPARSE_STOKES Apply prebuilt sparse pair maps.
%
% This is the solve-grid, cmap/self-correct mobility peanut transform used
% by matvec_mob_peanut_big_sparse. Fine/image source outputs are empty in
% this v1 path; postprocessing keeps using transform_mob_peanut_stokes.

if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('transform_mob_peanut_big_sparse_stokes:MissingBigSparse', ...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('transform_mob_peanut_big_sparse_stokes:UnsupportedGrid', ...
        'The big sparse transform is only supported on the solve grid.');
end

rvec_out = geom.rvec_out;
q = geom.q(:);
opt = geom.opt;
U = basis.U;
Y = basis.Y;
Lc = basis.Lc;

P = numel(q);
N_c = opt.N_c;
N_large = numel(rvec_out)/P;
PM = numel(rvec_out);
n_coarse = P*N_c;

lam_c_nonpx = zeros(n_coarse,1);
lam_c_nonpy = zeros(n_coarse,1);
lam_c_x = zeros(n_coarse,1);
lam_c_y = zeros(n_coarse,1);

% Keep the one-body source construction identical to the serial transform.
for i = 1:P
    coarse_ind = (i-1)*N_c+1:i*N_c;
    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    tau_mapped = Y{1}*(U{1}*[tau_particle_x; tau_particle_y]);

    lam_c_nonpx(coarse_ind) = tau_mapped(1:N_c);
    lam_c_nonpy(coarse_ind) = tau_mapped(N_c+1:end);

    tau_i_x = tau_mapped(1:N_c);
    tau_i_y = tau_mapped(N_c+1:end);
    tau_proj = [tau_i_x; tau_i_y] - Lc*[tau_i_x; tau_i_y];

    lam_c_x(coarse_ind) = tau_proj(1:N_c);
    lam_c_y(coarse_ind) = tau_proj(N_c+1:end);
end

lam_self_x = lam_c_x;
lam_self_y = lam_c_y;

lambda_self = [lam_self_x; lam_self_y];
use_direct_source_corr = logical(getOptField(opt, ...
    'big_sparse_direct_source_corr',false));
use_direct_u_corr = logical(getOptField(opt,'big_sparse_direct_u_corr',false));
use_structured_apply = logical(getOptField(opt, ...
    'big_sparse_structured_apply',false));
use_combined_pair_u = isfield(basis.big_sparse,'M_pair_nonp_u_corr');
need_pair_proj = ~use_direct_source_corr || ~use_direct_u_corr;
u_corr = [];
if need_pair_proj
    if use_combined_pair_u
        combined = basis.big_sparse.M_pair_nonp_u_corr*lambda_self;
        pair_rows = basis.big_sparse.pair_nonp_rows;
        pair_nonp = combined(1:pair_rows);
        u_corr = combined(pair_rows+1:end);
    elseif isfield(basis.big_sparse,'M_pair_nonp')
        pair_nonp = basis.big_sparse.M_pair_nonp*lambda_self;
    else
        error('transform_mob_peanut_big_sparse_stokes:MissingPairMaps', ...
            'Factored source or velocity correction requires M_pair_nonp.');
    end
    % Algebra per pair:
    %   pair_nonp = C_nonp*lambda_self
    %   pair_proj = P_pair*pair_nonp
    if use_structured_apply
        if ~isfield(basis.big_sparse,'P_pair')
            error('transform_mob_peanut_big_sparse_stokes:MissingProjector', ...
                'Structured apply requires basis.big_sparse.P_pair.');
        end
        pair_proj = applyDensePairProjector(pair_nonp, ...
            basis.big_sparse.P_pair);
    else
        if ~isfield(basis.big_sparse,'M_pair_projector')
            error('transform_mob_peanut_big_sparse_stokes:MissingProjector', ...
                'Sparse factored apply requires M_pair_projector.');
        end
        pair_proj = basis.big_sparse.M_pair_projector*pair_nonp;
    end
else
    pair_proj = [];
    pair_nonp = [];
end

if use_direct_source_corr
    % Diagnostic path: M_source_corr forms P_pair*C_nonp explicitly. The
    % default factored source correction is more stable for large sources.
    if ~isfield(basis.big_sparse,'M_source_corr')
        error('transform_mob_peanut_big_sparse_stokes:MissingDirectSourceCorr', ...
            'opt.big_sparse_direct_source_corr=1 requires M_source_corr.');
    end
    corr = basis.big_sparse.M_source_corr*lambda_self;
else
    if use_structured_apply
        if ~isfield(basis.big_sparse,'source_scatter_rows')
            error('transform_mob_peanut_big_sparse_stokes:MissingSourceRows', ...
                'Structured source correction requires source_scatter_rows.');
        end
        corr = scatterSourceCorrections(basis.big_sparse.source_scatter_rows, ...
            pair_proj,pair_nonp,4*N_c,n_coarse);
    else
        if ~isfield(basis.big_sparse,'M_source_scatter')
            error('transform_mob_peanut_big_sparse_stokes:MissingSourceScatter', ...
                'Factored source correction requires M_source_scatter.');
        end
        corr = basis.big_sparse.M_source_scatter*[pair_proj; pair_nonp];
    end
end

lam_c_x = lam_self_x + corr(1:n_coarse);
lam_c_y = lam_self_y + corr(n_coarse+1:2*n_coarse);
lam_c_nonpx = lam_c_nonpx + corr(2*n_coarse+1:3*n_coarse);
lam_c_nonpy = lam_c_nonpy + corr(3*n_coarse+1:4*n_coarse);

if use_direct_u_corr
    % Diagnostic path: M_u_corr forms Ucross - Ecolloc*C_proj explicitly,
    % which can lose precision compared with the default factored apply.
    if ~isempty(u_corr)
        % Computed together with pair_nonp by M_pair_nonp_u_corr.
    elseif ~isfield(basis.big_sparse,'M_u_corr')
        error('transform_mob_peanut_big_sparse_stokes:MissingDirectUCorr', ...
            'opt.big_sparse_direct_u_corr=1 requires M_u_corr.');
    else
        u_corr = basis.big_sparse.M_u_corr*lambda_self;
    end
elseif isfield(basis.big_sparse,'M_u_cross') && ...
        isfield(basis.big_sparse,'M_u_peanut')
    % Default path: u_corr = Ucross*lambda_self - Ecolloc*pair_proj.
    u_corr = basis.big_sparse.M_u_cross*lambda_self - ...
        basis.big_sparse.M_u_peanut*pair_proj;
elseif isfield(basis.big_sparse,'M_u_corr')
    u_corr = basis.big_sparse.M_u_corr*lambda_self;
else
    error('transform_mob_peanut_big_sparse_stokes:MissingUCorrMap', ...
        'Missing both factored u-correction maps and direct M_u_corr.');
end

lam_f_x = [];
lam_f_y = [];
rimage_k = cell(P,1);
end

function pair_proj = applyDensePairProjector(pair_nonp,P_pair)
if isempty(pair_nonp)
    pair_proj = pair_nonp;
    return
end
pair_rows = size(P_pair,1);
n_pairs = numel(pair_nonp)/pair_rows;
pair_proj = P_pair*reshape(pair_nonp,pair_rows,n_pairs);
pair_proj = pair_proj(:);
end

function corr = scatterSourceCorrections(rows,pair_proj,pair_nonp, ...
    pair_rows,n_coarse)
n_source_rows = 4*n_coarse;
if isempty(rows)
    corr = zeros(n_source_rows,1);
    return
end
n_pairs = numel(pair_nonp)/pair_rows;
pair_proj = reshape(pair_proj,pair_rows,n_pairs);
pair_nonp = reshape(pair_nonp,pair_rows,n_pairs);
vals = reshape([pair_proj; pair_nonp],[],1);
corr = accumarray(rows,vals,[n_source_rows,1]);
end
