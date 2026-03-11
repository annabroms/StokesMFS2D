function res = applyPmat_peanut(vel,rvec_in,rvec_out,Sinv,q,Zi,Yi,rbase_in_c,...
    rbase_in_f,rbase_out_f,refine,rimage,nimage,opt,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap)
%APPLYPMAT  Apply long-range preconditioning projection matrix to a velocity field.
%
%   res = APPLYPMAT_PEANUT(vel, rvec_in, rvec_out, Rinv, Zi, Yi, opt)
%
%   Applies the projector 
%
%      %       P = I - A * (Z * (S \ I) * Y')
%
%   to the input velocity vector, with Z = blkdiag(Zi) and Y = blkdiag(Yi).
%   This corresponds to projecting out contributions from the coarse subspace 
%   used in long-range preconditioning applied to A, mapping mu to boundary
%   velocities (a square matrix constructed via two-body preconditioning with peanut compression).
%   
%
%   INPUTS:
%     vel       - 2PM×1 velocity vector at the M boundary points of P bodies.
%     rvec_in   - PN complex valued array of source points for the projection flow.
%     rvec_out  - PM complex valued array of target points where the projection flow is evaluated.
%     Zi         - (2M×k) matrix mapping coarse coefficients to proxy forces for a single body. k is
%                 set by the number of coarse basis functions per body
%     Yi         - (2M×k) matrix such that Y' maps surface flow to coarse space for a single body.
%     Sinv      - (Pk×Pk) inverse of coarse interaction matrix S.  
%     opt       - Struct with fmm and kmax flags etc
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
%   See also: applyQmat, getCoarseSource, get_long_range_precond,
%   getVelocityField
%
% Anna Broms Oct 19, 2025

mu = getCoarseMu(vel,Sinv,Zi,Yi,opt.db,opt.P,opt.N_c,opt.a_c);

%compute velocities using these source strengths

%remember that everything so far is implemented without images in mind
debug = 0; 

proj_vel = matvec_2B_peanut(mu,rbase_in_c,rbase_in_f,rvec_in,rbase_out_f,...
    refine,rimage,nimage,opt,rvec_out,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,debug);
%UU{1} = UU{1}'; %stupid
%pair_points = [zeros(opt.P,1) opt.N_c*opt.a_c*ones(opt.P,1)];
%proj_vel = matvec_2D_Stokes(mu,rvec_in,rvec_out,[],[],q,UU,YY,pair_points,opt.s);

res = vel-proj_vel;


end