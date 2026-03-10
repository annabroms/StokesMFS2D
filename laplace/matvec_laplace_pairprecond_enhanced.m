function res = matvec_laplace_pairprecond_enhanced(tau,geom,basis)
%MATVEC_LAPLACE_PAIRPRECOND_ENHANCED Matrix-vector action for Laplace 2B preconditioner.
%
% Syntax:
%   res = matvec_laplace_pairprecond_enhanced(tau,geom,basis)
%
% Anna Broms, Mar 2026

rbase_in_c = geom.rbase_in_c;
rvec_out = geom.rvec_out;
q = geom.q;
opt = geom.opt;

P = length(q);
N_large = length(rvec_out)/P;
N_c = opt.N_c;

[rvec_in,coarse_ind,lam_all,u_corr,lam_c] = getPairTransformationLaplace(tau,geom,basis);

res = laplaceSingleLayerField(rvec_in,rvec_out,lam_all,opt.use_fmm);

rout = rvec_out(1:N_large)-q(1);
Nii = lapSLPmat(rbase_in_c,rout);

for i = 1:P
    idx = coarse_ind{i};
    uii = Nii*lam_c(idx);
    block = (i-1)*N_large+1:i*N_large;
    res(block) = res(block)-uii;
end

res = res-u_corr;
res = res+tau;

if length(res) ~= P*ceil(opt.a_c*N_c)
    % No-op guard; keeps diagnostics close to Stokes style.
end

end
