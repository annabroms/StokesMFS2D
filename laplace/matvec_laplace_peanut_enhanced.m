function res = matvec_laplace_peanut_enhanced(tau,geom,basis)
%MATVEC_LAPLACE_PEANUT_ENHANCED Matrix-vector action for Laplace peanut preconditioner.
%
% Syntax:
%   res = matvec_laplace_peanut_enhanced(tau,geom,basis)

rbase_in_c = geom.rbase_in_c;
rvec_in = geom.rvec_in;
rvec_out = geom.rvec_out;
opt = geom.opt;
q = geom.q;

P = length(q);
N_large = length(rvec_out)/P;

[lam_c,lam_self,~,~,u_corr] = transform_laplace_peanut(tau,geom,basis);

res = laplaceSingleLayerField(rvec_in,rvec_out,lam_c,opt.use_fmm);
res = res + u_corr;

rout = rvec_out(1:N_large)-q(1);
Nii = lapSLPmat(rbase_in_c,rout);

for i = 1:P
    idx = (i-1)*opt.N_c+1:i*opt.N_c;
    block = (i-1)*N_large+1:i*N_large;
    uii = Nii*lam_self(idx);
    res(block) = res(block)-uii;
end

res = res + tau;

end
