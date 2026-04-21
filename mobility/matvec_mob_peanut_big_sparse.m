function res = matvec_mob_peanut_big_sparse(tau,geom,basis)
%MATVEC_MOB_PEANUT_BIG_SPARSE Mobility peanut matvec using sparse pair maps.
%
% This v1 matvec is only for the GMRES solve grid with opt.cmap=1,
% opt.self_correct=1, opt.use_dense=1, and opt.use_big_sparse=1.

rbase_in_c = geom.rbase_in_c;
rvec_in = geom.rvec_in;
opt = geom.opt;
rvec_out = geom.rvec_out;
rcheck = geom.rcheck;
q = geom.q(:);
use_fmm = logical(getOptField(opt,'use_fmm',true));

validateBigSparseMatvecOptions(geom,basis);

P = numel(q);
PM = numel(rvec_out);
PM_check = numel(rcheck);
N_large = PM_check/P;
N_c = opt.N_c;
mu = 1;

[lam_c_x, ~, lam_self_x, ~, lam_c_y, ~, lam_self_y, ~, u_corr] = ...
    transform_mob_peanut_big_sparse_stokes(tau,geom,basis);

res = zeros(2*PM_check,1);

% Correct identity blocks on the solve grid, as in matvec_mob_peanut_enhanced.
rout = rvec_out(1:N_large)-q(1);
Nii = stokSLPmat(rbase_in_c,rout,mu);
for i = 1:P
    coarse_ind = (i-1)*N_c+1:i*N_c;
    tau_xy = [lam_self_x(coarse_ind); lam_self_y(coarse_ind)];
    uii = Nii*tau_xy;

    res((i-1)*N_large+1:i*N_large) = ...
        res((i-1)*N_large+1:i*N_large) - uii(1:end/2);
    res((i-1)*N_large+PM+1:i*N_large+PM) = ...
        res((i-1)*N_large+PM+1:i*N_large+PM) - uii(end/2+1:end);
end

res = res + tau;

vel = getVelocityField(rvec_in,rcheck,lam_c_x,lam_c_y,use_fmm);
res = res + vel + u_corr;
end

function validateBigSparseMatvecOptions(geom,basis)
opt = geom.opt;
if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('matvec_mob_peanut_big_sparse:MissingBigSparse', ...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('matvec_mob_peanut_big_sparse:UnsupportedGrid', ...
        'The big sparse mobility matvec is only supported on the solve grid.');
end
if ~logical(getOptField(opt,'cmap',false))
    error('matvec_mob_peanut_big_sparse:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.cmap=1.');
end
if ~logical(getOptField(opt,'self_correct',false))
    error('matvec_mob_peanut_big_sparse:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.self_correct=1.');
end
if ~logical(getOptField(opt,'use_dense',false))
    error('matvec_mob_peanut_big_sparse:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.use_dense=1.');
end
end
