function res = matvec_2D_peanut_enhanced(tau,geom,basis)
%MATVEC_2D_PEANUT_ENHANCED Matrix-vector action for resistance peanut preconditioner based on Stokeslet souces only.
%
% Syntax:
%   res = matvec_2D_peanut_enhanced(tau,geom,basis)
%
% Inputs:
%   tau   - Stacked boundary data [tau_x; tau_y] on outer collocation points.
%   geom  - Geometry/problem struct with fields:
%           rbase_in_c, rvec_in, opt, rvec_out, rcheck, q, pairs,
%           rbase_in_f, refine, rimage_vec.
%   basis - Precomputed basis struct with fields:
%           U, Y, Upf, Ypf, DC_all, YC_all, Cmap.
%
% Output:
%   res   - Matvec result at collocation/check points [u_x; u_y].

rbase_in_c = geom.rbase_in_c;
rvec_in = geom.rvec_in;
opt = geom.opt;
rvec_out = geom.rvec_out;
q = geom.q;

P = length(q);
N_large = length(rvec_out)/P;
PM = length(rvec_out);
N_c = opt.N_c;
mu = 1;

if isfield(opt,'use_fmm_velocity')
    use_fmm_velocity = opt.use_fmm_velocity;
else
    use_fmm_velocity = true;
end

% Transform coarse boundary data to compressed coarse/fine source strengths.
[lam_stokes_x, lam_self_x, ~, lam_stokes_y, lam_self_y, ~, u_corr] = ...
    transform_peanut_stokes(tau,geom,basis);

% Evaluate flow induced by compressed source strengths.
res = getVelocityField(rvec_in,rvec_out,lam_stokes_x,lam_stokes_y,use_fmm_velocity);

% Replace compressed pair blocks by local fine representation.
res = res + u_corr;

% Subtract duplicated self-contribution and enforce identity on diagonal.
rout = rvec_out(1:N_large)-q(1);
Nii = singleLayer(rbase_in_c,rout,mu);

for i = 1:P
    coarse_ind = (i-1)*N_c+1:i*N_c;
    tau_xy = [lam_self_x(coarse_ind); lam_self_y(coarse_ind)];
    uii = Nii*tau_xy;

    % x-component block
    res((i-1)*N_large+1:i*N_large) = ...
        res((i-1)*N_large+1:i*N_large) - uii(1:end/2) + tau((i-1)*N_large+1:i*N_large);

    % y-component block
    res((i-1)*N_large+PM+1:i*N_large+PM) = ...
        res((i-1)*N_large+PM+1:i*N_large+PM) - uii(end/2+1:end) + ...
        tau((i-1)*N_large+PM+1:i*N_large+PM);
end

end
