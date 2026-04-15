function res = matvec_res_3B_enhanced(tau,geom,basis)
%MATVEC_RES_3B_ENHANCED Matrix-vector action for resistance triplet preconditioner.

rbase_in_c = geom.rbase_in_c;
opt = geom.opt;
rvec_out = geom.rvec_out;
q = geom.q;

[rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y,u_corr] = ...
    getTripletTransformationStokes(tau,geom,basis);

use_fmm = logical(getOptField(opt,'use_fmm',true));
res = getVelocityField(rvec_in,rvec_out,tau_stokes_x,tau_stokes_y,use_fmm);

P = opt.P;
PM = numel(rvec_out);
N_large = PM/P;
Nii = stokSLPmat(rbase_in_c,rvec_out(1:N_large)-q(1),1);

for i = 1:P
    tau_xy = [tau_stokes_x(coarse_ind{i}); tau_stokes_y(coarse_ind{i})];
    uii = Nii*tau_xy;

    segx = (i-1)*N_large+1:i*N_large;
    segy = PM + segx;
    res(segx) = res(segx) - uii(1:N_large);
    res(segy) = res(segy) - uii(N_large+1:end);
end

res = res - u_corr + tau;

end
