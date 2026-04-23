function [lam_c_x, lam_c_y, lambda_self_blocks, u_corr, ...
    lam_c_nonpx, lam_c_nonpy] = ...
    transform_mob_peanut_big_sparse_stokes(tau,geom,basis)
%TRANSFORM_MOB_PEANUT_BIG_SPARSE_STOKES Apply prebuilt sparse pair maps.
%
% This is the fast solve-grid, cmap/self-correct mobility peanut transform
% used by matvec_peanut_big_sparse. The first four outputs match the
% resistance big-sparse transform: coarse x/y sources, projected one-body
% sources in block form, and the close-pair velocity correction.
% Boundary postprocessing can additionally request the nonprojected coarse
% one-body sources.

if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('transform_mob_peanut_big_sparse_stokes:MissingBigSparse', ...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('transform_mob_peanut_big_sparse_stokes:UnsupportedGrid', ...
        'The big sparse transform is only supported on the solve grid.');
end

rvec_out = geom.rvec_out;
opt = geom.opt;
U = basis.U;
Y = basis.Y;
Lc = basis.Lc;

P = getOptField(opt,'P',numel(geom.q));
N_c = opt.N_c;
PM = numel(rvec_out);
N_large = PM/P;
n_coarse = P*N_c;

% 1-body: recover sources from boundary unknowns
% Keep the one-body source construction in the same U{1} then Y{1} order
% as the serial transform, but apply it to all particle columns at once.
tau_x = reshape(tau(1:PM),N_large,P);
tau_y = reshape(tau(PM+1:2*PM),N_large,P);
tau_blocks = [tau_x; tau_y];
tau_mapped_blocks = Y{1}*(U{1}*tau_blocks);

% project
tau_proj_blocks = tau_mapped_blocks - Lc*tau_mapped_blocks;
lambda_self_blocks = tau_proj_blocks;

lam_c_nonpx = reshape(tau_mapped_blocks(1:N_c,:),[],1);
lam_c_nonpy = reshape(tau_mapped_blocks(N_c+1:2*N_c,:),[],1);
lam_c_x = reshape(tau_proj_blocks(1:N_c,:),[],1);
lam_c_y = reshape(tau_proj_blocks(N_c+1:2*N_c,:),[],1);

lambda_self = [reshape(lambda_self_blocks(1:N_c,:),[],1); ...
    reshape(lambda_self_blocks(N_c+1:2*N_c,:),[],1)];
if ~isfield(basis.big_sparse,'M_pair_nonp')
    error('transform_mob_peanut_big_sparse_stokes:MissingPairMap', ...
        'The big sparse transform requires M_pair_nonp.');
end
if ~isfield(basis.big_sparse,'P_pair')
    error('transform_mob_peanut_big_sparse_stokes:MissingProjector', ...
        'The big sparse transform requires basis.big_sparse.P_pair.');
end
if ~isfield(basis.big_sparse,'source_scatter_rows')
    error('transform_mob_peanut_big_sparse_stokes:MissingSourceRows', ...
        'The big sparse transform requires source_scatter_rows.');
end

% Fixed source-correction algebra:
%   pair_nonp = C_nonp*lambda_self
%   pair_proj = P_pair*pair_nonp
%   corr = scatter([pair_proj; pair_nonp])
pair_nonp = basis.big_sparse.M_pair_nonp*lambda_self;
pair_proj = applyDensePairProjector(pair_nonp,basis.big_sparse.P_pair);
corr = scatterSourceCorrections(basis.big_sparse.source_scatter_rows, ...
    pair_proj,pair_nonp,4*N_c,n_coarse);

lam_c_x = lam_c_x + corr(1:n_coarse);
lam_c_y = lam_c_y + corr(n_coarse+1:2*n_coarse);
% Keep lam_c_nonp on the same convention as transform_mob_peanut_stokes:
% one-body nonprojected sources only. Pair rigid-motion contributions are
% recovered separately from Cmap_FU/M_rbm_corr in postprocessing.

use_direct_u_corr = logical(getOptField(opt,'big_sparse_direct_u_corr',true));
if use_direct_u_corr
    if ~isfield(basis.big_sparse,'M_u_corr')
        error('transform_mob_peanut_big_sparse_stokes:MissingUCorr', ...
            'opt.big_sparse_direct_u_corr=1 requires M_u_corr.');
    end
    u_corr = basis.big_sparse.M_u_corr*lambda_self;
else
    if ~isfield(basis.big_sparse,'M_u_cross') || ...
            ~isfield(basis.big_sparse,'M_u_peanut')
        error('transform_mob_peanut_big_sparse_stokes:MissingFactoredUCorr', ...
            ['opt.big_sparse_direct_u_corr=0 requires M_u_cross ', ...
             'and M_u_peanut.']);
    end
    u_corr = basis.big_sparse.M_u_cross*lambda_self - ...
        basis.big_sparse.M_u_peanut*pair_proj;
end


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
