function res = matvec_lap_peanut_big_sparse(tau,geom,basis)
%MATVEC_LAP_PEANUT_BIG_SPARSE Laplace peanut matvec using sparse pair maps.

validateBigSparseMatvecOptions(geom,basis);

rbase_in_c = geom.rbase_in_c;
rvec_in = geom.rvec_in;
rvec_out = geom.rvec_out;
opt = geom.opt;
q = geom.q(:);
pairs = geom.pairs;
rcheck = geom.rcheck;

P = numel(q);
N_large = numel(rvec_out)/P;

[lam_c,lam_self,~,~,u_corr,pair_qv_nonp] = ...
    transform_lap_peanut_big_sparse(tau,geom,basis);

res = lapSLPfield(rvec_in,rcheck,lam_c,opt.use_fmm);
res = res + u_corr;

if isfield(opt,'project_charge') && logical(opt.project_charge)
    has_neigh = sort(unique(pairs(:)));
    for ii = 1:numel(has_neigh)
        k = has_neigh(ii);
        block = (k-1)*N_large+1:k*N_large;
        res(block) = res(block) + pair_qv_nonp(k);
    end
end

if isfield(basis,'Nii') && ~isempty(basis.Nii)
    Nii = basis.Nii;
else
    rout = rvec_out(1:N_large)-q(1);
    Nii = lapSLPmat(rbase_in_c,rout);
end

for i = 1:P
    idx = (i-1)*opt.N_c+1:i*opt.N_c;
    block = (i-1)*N_large+1:i*N_large;
    uii = Nii*lam_self(idx);
    res(block) = res(block) - uii;
end

res = res + tau;
end

function validateBigSparseMatvecOptions(geom,basis)
opt = geom.opt;
if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('matvec_lap_peanut_big_sparse:MissingBigSparse', ...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('matvec_lap_peanut_big_sparse:UnsupportedGrid', ...
        'The Laplace big-sparse matvec is only supported on the solve grid.');
end
if ~logical(getOptField(opt,'cmap',false))
    error('matvec_lap_peanut_big_sparse:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.cmap=1.');
end
end
