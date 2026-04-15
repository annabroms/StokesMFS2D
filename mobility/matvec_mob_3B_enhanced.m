function res = matvec_mob_3B_enhanced(tau,geom,basis)
%MATVEC_MOB_3B_ENHANCED Matrix-vector action for mobility triplet preconditioner.

rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
opt = geom.opt;
rvec_out = geom.rvec_out;
rcheck = geom.rcheck;
q = geom.q;
triplets = geom.triplets;

is_collocation = isequal(rcheck,rvec_out);
[rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y, ...
    tau_stokes_nonpx,tau_stokes_nonpy,tau_stokes_e_nonpx,tau_stokes_e_nonpy, ...
    rimage_k] = getMobTripletTransformationStokes(tau,geom,basis);

use_fmm = logical(getOptField(opt,'use_fmm',true));
res = getVelocityField(rvec_in,rcheck,tau_stokes_x,tau_stokes_y,use_fmm);
if ~is_collocation
    return
end

P = opt.P;
PM = numel(rvec_out);
N_large = PM/P;
N_c = opt.N_c;
N_f = opt.N_f;
rbase_out_rel = rvec_out(1:N_large)-q(1);

has_trip = sort(unique(triplets(:)));
offset = P*N_c;
for it = 1:numel(has_trip)
    pid = has_trip(it);

    fine_idx = offset + (1:N_f);
    offset = offset + N_f;
    fcx = tau_stokes_nonpx(fine_idx);
    fcy = tau_stokes_nonpy(fine_idx);
    bcvec_f = applyBKt2D(rbase_out_rel,0,rbase_in_f,0,fcx,fcy);

    if isempty(rimage_k{pid})
        bcvec_e = zeros(2*N_large,1);
    else
        bcvec_e = applyBKt2D(rbase_out_rel,0,rimage_k{pid},q(pid), ...
            tau_stokes_e_nonpx{pid},tau_stokes_e_nonpy{pid});
    end

    segx = (pid-1)*N_large+1:pid*N_large;
    segy = PM + segx;
    res(segx) = res(segx) + bcvec_f(1:N_large) + bcvec_e(1:N_large);
    res(segy) = res(segy) + bcvec_f(N_large+1:end) + bcvec_e(N_large+1:end);
end

if offset ~= numel(tau_stokes_nonpx)
    error('matvec_mob_3B_enhanced:IndexMismatch', ...
        'Unexpected fine-source indexing mismatch.');
end

Nii = stokSLPmat(rbase_in_c,rbase_out_rel,1);
for i = 1:P
    tau_xy = [tau_stokes_x(coarse_ind{i}); tau_stokes_y(coarse_ind{i})];
    uii = Nii*tau_xy;

    segx = (i-1)*N_large+1:i*N_large;
    segy = PM + segx;
    res(segx) = res(segx) - uii(1:N_large) + tau(segx);
    res(segy) = res(segy) - uii(N_large+1:end) + tau(segy);
end

end
