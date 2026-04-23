function res = matvec_peanut_big_sparse(tau,geom,basis)
%MATVEC_PEANUT_BIG_SPARSE Stokes peanut matvec using sparse pair
%maps. Used both for resistance and mobility.
% This is an optimized version of the matvecs that use a loop over near
% contacts in matvec_res_peanut_enhanced and matvec_mob_peanut_enhanced. 
% It also uses batched versions of  the 1-body pseudoinverse and identity correction.
% For this, see comparison against the loop over all bodies in
% apr23_resistance_one_body_batch_compare.m

validateBigSparseMatvecOptions(geom,basis);

rvec_in = geom.rvec_in;
rvec_out = geom.rvec_out;
opt = geom.opt;
use_fmm = opt.use_fmm;
P = numel(geom.q);
PM = numel(rvec_out);
N_large = PM/P;

% Recover compressed source strengths and the close-pair velocity fix.
if opt.resistance
    [lam_c_x, lam_c_y, lambda_self_blocks, u_corr] = ...
        transform_res_peanut_big_sparse_stokes(tau,geom,basis);
else
    [lam_c_x, lam_c_y, lambda_self_blocks, u_corr] = ...
        transform_mob_peanut_big_sparse_stokes(tau,geom,basis);
end

% Evaluate the global coarse Stokeslet field 
res = getVelocityField(rvec_in,rvec_out,lam_c_x,lam_c_y,use_fmm);
% add pair velocity corrections to enforce diagonal blocks
res = res + u_corr;

% Replace each diagonal one-body block by the identity contribution.

uii_blocks = basis.Nii * lambda_self_blocks;

res(1:PM) = res(1:PM) - reshape(uii_blocks(1:N_large,:),[],1);
res(PM+1:2*PM) = res(PM+1:2*PM) - ...
    reshape(uii_blocks(N_large+1:2*N_large,:),[],1);

res = res + tau;

end

function validateBigSparseMatvecOptions(geom,basis)
opt = geom.opt;
if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('matvec_peanut_big_sparse:MissingBigSparse', ...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('matvec_peanut_big_sparse:UnsupportedGrid', ...
        'The resistance big sparse matvec is only supported on the solve grid.');
end
if ~logical(getOptField(opt,'cmap',false))
    error('matvec_peanut_big_sparse:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.cmap=1.');
end
if ~logical(getOptField(opt,'self_correct',false))
    error('matvec_peanut_big_sparse:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.self_correct=1.');
end

if ~logical(getOptField(opt,'use_dense',false))
    error('matvec_peanut_big_sparse:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.use_dense=1.');
end
end
