function res = applyPmat_mu(vel, rvec_in, rvec_out, Sinv, q, Zi, Yi, Uii, Yii, pair_points, opt)
%APPLYPMAT_MU Apply long-range preconditioning projection matrix to a velocity field.
%
%   res = APPLYPMAT_MU(vel, rvec_in, rvec_out, Sinv, q, Zi, Yi, ...
%                      rimage, nimage, Uii, Yii, pair_points, opt)
%
%   Applies the projection operator
%
%       P = I - A * (Z * (S^{-1}) * Y')
%
%   to a given velocity field, removing the coarse-space components
%   associated with the long-range preconditioner in the MFS mobility solver.
%   The operator A corresponds to the boundary-to-boundary mapping obtained with 
%   one-body preconditioning from the right of the original least squares problem,
%   while Z and Y define the coarse subspace bases per particle.
%
%   INPUTS:
%       vel         - 2PM × 1 vector of boundary velocities for P particles, 
%                     with M discretization points per boundary.
%       rvec_in     - Complex-valued PN × 1 array of source (proxy) points for projection flow.
%       rvec_out    - Complex-valued PM × 1 array of target (collocation) points.
%       Sinv        - (Pk × Pk) inverse of the coarse interaction matrix S 
%                     
%       q           - Complex valued array of the P particle center coordinates.
%       Zi          - (2M × k) matrix mapping coarse coefficients to proxy sources 
%                     for one particle (defines the Z block).
%       Yi          - (2M × k) matrix whose transpose maps surface velocity to 
%                     coarse coefficients for one particle (defines the Y block).
%       Uii         - self-interaction matrix (one-body operator for each
%                     particle) that together with Yii builds the pseudoinverse
%       Yii         
%       pair_points - matric containing the number of image points and collocation points for body k on row k s
%       opt         - Struct with options controlling:
%                       • opt.use_fmm  – logical flag for FMM acceleration
%                       • opt.images   – include image interactions if true
%
%   OUTPUT:
%     res       - Projected velocity vector res = P vel
%                 where lambda are the coarse source strengths.
%
%   NOTES:
%     - The projector removes the coarse flow contribution defined by the span
%       of Z, and used in the long-range preconditioner.
%     - The matrices Z and Y are sorted as [x x ... x y y ... y ]
%     - Currently implemented without image system contributions.
%
%   See also: applyQmat_mu, getCoarseMu, get_long_range_precond,
%   getVelocityField, 
%
% Anna Broms Oct 19, 2025

mu = getCoarseMu(vel,Sinv,Zi,Yi,opt.db,opt.P,opt.N_c,opt.a_c);

%compute velocities using these source strengths

%remember that everything so far is implemented without images in mind

proj_vel = matvec_2D_Stokes(mu,rvec_in,rvec_out,[],[],q,Uii,Yii,pair_points,opt.s);

res = vel-proj_vel;


end