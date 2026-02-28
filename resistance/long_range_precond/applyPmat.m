function res = applyPmat(vel,rvec_in,rvec_out,Sinv,q,Ny,Mx,Zi,Yi,opt)
%APPLYPMAT  Apply long-range preconditioning projection matrix to a velocity field.
%
%   res = APPLYPMAT(vel, rvec_in, rvec_out, Rinv, Zi, Yi, opt)
%
%   Applies the projector 
%
%       P = I - G * (Z * (R \ I) * Y')
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

lambda = getCoarseSource(vel,Sinv,Zi,Yi,opt.db,opt.P,opt.N_c,opt.a_c);

%compute velocities using these source strengths

%remember that everything so far is implemented without images in mind

%mask out elements here
if opt.mask
    Nc = opt.N_c;
    a = opt.a_c;
    P = opt.P;
    vel_k = zeros(P*Nc*a*2,1);
    proj_vel = zeros(P*Nc*a*2,1);
    for k = 1:P
        vel_k(:) = 0;
        lambda_x = lambda(Nc*(k-1)+1:Nc*k);
        lambda_y = lambda(Nc*(k-1)+P*Nc+1:Nc*k+P*Nc);
        proj_k = getVelocityField(rvec_in(Nc*(k-1)+1:Nc*k),rvec_out,lambda_x,lambda_y,[],[],[],[],[],[], []);
        d = abs(q-q(k));
        ind = find(d>(2+opt.cut_off));
        ind_diff = setdiff(1:P,ind);
        keep_ind = []; 
        for l = 1:length(ind_diff)
            i = ind_diff(l); 
            keep_ind = [keep_ind; ((i-1)*Nc*a+1:i*Nc*a)'];
        end
        vel_k(keep_ind) = proj_k(keep_ind);
        vel_k(keep_ind+P*Nc*a) = proj_k(keep_ind+P*Nc*a);
    
        proj_vel = proj_vel+vel_k;
    
    
    end
else
    proj_vel = getVelocityField(rvec_in,rvec_out,lambda(1:end/2),lambda(end/2+1:end),[],[],[],[],[],[], []);
end

res = vel-proj_vel;


end