function res = applyPmat_mob(vel,rvec_in,rvec_out,L,Lr,Sinv,Zi,Yi,opt)
%APPLYPMAT_MOB  Apply long-range preconditioning projection matrix to a velocity field.
%
%   res = APPLYPMAT_MOB(vel, rvec_in, rvec_out, Rinv, Zi, Yi, opt)
%
%   Applies the projector 
%
%       P = I - G * (Z * (S \ I) * Y')
%
%   to the input velocity vector, with Z = blkdiag(Zi) and Y = blkdiag(Yi).
%   This corresponds to projecting out contributions from the coarse subspace 
%   used in long-range preconditioning.
%   
%
%   INPUTS:
%     vel       - 2PM×1 velocity vector at the M boundary points of P bodies.
%     rvec_in   - PN complex valued array of source points for the projection flow.
%     rvec_out  - PM complex valued array of target points where the projection flow is evaluated.
%     Zi         - (2N×k) matrix mapping coarse coefficients to proxy forces for a single body. k is
%                 set by the number of coarse basis functions per body
%     Yi         - (2M×k) matrix such that Y' maps surface flow to coarse space for a single body.
%     Sinv      - (Pk×Pk) inverse of coarse interaction matrix R.  
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
%   See also: applyQmat_mob, getCoarseSource, get_long_range_precond_mob,
%   getVelocityField
%
% Anna Broms Oct 19, 2025

lambda = getCoarseSource(vel,Sinv,Zi,Yi,opt.db,opt.P,opt.N_c,opt.a_c);

%compute velocities using these source strengths

P = opt.P;
Nc = opt.N_c; 
a = opt.a_c;

proj_L = zeros(P*Nc*2,1);
Lr_vel = zeros(P*Nc*a*2,1);
for k = 1:opt.P
    vel_k = [lambda(Nc*(k-1)+1:Nc*k); lambda(P*Nc+Nc*(k-1)+1:Nc*k+P*Nc)];
    proj_vel_k = vel_k-L*vel_k; %project by I-L
    proj_L(Nc*(k-1)+1:Nc*k) = proj_vel_k(1:end/2);
    proj_L(Nc*(k-1)+P*Nc+1:Nc*k+P*Nc) = proj_vel_k(end/2+1:end);

    Lr_vel_k = Lr*vel_k; % Lr lambda
    Lr_vel(Nc*a*(k-1)+1:Nc*a*k) = Lr_vel_k(1:end/2);
    Lr_vel(Nc*a*(k-1)+P*Nc*a+1:Nc*a*k+P*Nc*a) = Lr_vel_k(end/2+1:end);


end

%remember that everything so far is implemented without images in mind
proj_vel = getVelocityField(rvec_in,rvec_out,proj_L(1:end/2),proj_L(end/2+1:end),[],[],[],[],[],[], []);



res = vel-proj_vel-Lr_vel;


end