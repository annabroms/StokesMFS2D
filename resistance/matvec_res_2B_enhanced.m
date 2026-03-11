function res = matvec_res_2B_enhanced(tau,geom,basis)
%MATVEC_RES_2B_ENHANCED Matrix-vector action for resistance 2-body preconditioned solve.
%
% Syntax:
%   res = matvec_res_2B_enhanced(tau,geom,basis)
%
% Inputs:
%   tau   - Stacked boundary data [tau_x; tau_y] on outer collocation points.
%   geom  - Geometry/problem struct with fields:
%           rbase_in_c, rbase_in_f, refine, rimage_vec, opt, rvec_out, q, pairs.
%   basis - Precomputed basis struct with fields:
%           U, Y, Upf, Ypf.
%
% Output:
%   res   - Matvec result at collocation points [u_x; u_y].

% See also: solve_res_2B_enhanced, getPairTransformationStokes.

rbase_in_c = geom.rbase_in_c;
opt = geom.opt;
rvec_out = geom.rvec_out;
q = geom.q;

% Transform data at collocation points to source strengths.
[rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y,u_corr] = ...
    getPairTransformationStokes(tau,geom,basis);

P = opt.P;
N_large = length(rvec_out)/P;
PM = length(rvec_out);
mu = 1;

%% Velocity field from all Stokeslet sources
res = getVelocityField(rvec_in,rvec_out,tau_stokes_x,tau_stokes_y);

two_corr = true; % identity correction for pair blocks

% Subtract self-interactions and replace diagonal with identity.
rout = rvec_out(1:N_large)-q(1,:);
Nii = singleLayer(rbase_in_c,rout,mu);

for i = 1:P
    tau_xy = [tau_stokes_x(coarse_ind{i}); tau_stokes_y(coarse_ind{i})];
    uii = Nii*tau_xy;

    % subtract contribution in x
    res((i-1)*N_large+1:i*N_large) = ...
        res((i-1)*N_large+1:i*N_large)-uii(1:end/2);

    % subtract contribution in y
    res((i-1)*N_large+1+PM:i*N_large+PM) = ...
        res((i-1)*N_large+1+PM:i*N_large+PM)-uii(end/2+1:end);
end

if two_corr
    res = res-u_corr;
end

res = res+tau;

end
