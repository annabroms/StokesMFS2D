function lambda_fine = applyQmat(vel,rvec_in,rvec_out,Rinv,Nx,Ny,Mx,Zi,Yi,opt)
%APPLYQMAT  Apply long-range preconditioning projection matrix on sources
%
%   lambda_fine = APPLYQMAT(vel, rvec_in, rvec_out, Rinv, Zi, Yi, opt)
%
%   Applies the projection operator 
%%   
%       Q = I - (Z * (R \ I) * Y') * G
%
%   to the solution vector lambda, with Z = blkdiag(Zi) and Y = blkdiag(Yi)
%
%   INPUTS:
%     lambda    - 3PN×1 vector of proxy source strengths 
%     rvec_in   - PN complex valued array of source points for evaluating the proxy flow.
%     rvec_out  - PM complex valued array of target surface points where the flow is evaluated.
%     Rinv      - (Pk×Pk) inverse of the coarse interaction matrix R,
%                 where k is the number of coarse basis functions per body
%     Y         - (2PM×Pk) matrix mapping from surface flow to coarse space.
%     Z         - (2PN×Pk) matrix mapping from proxy source space to coarse space.
%     opt       - Struct with flow options (e.g., FMM flags, kernel type).
%
%   OUTPUT:
%     lamda_fine       - Projected vector lambda_fine = Q * lambda, with coarse source contributions removed.
%
%   NOTES:
%     - Unlike APPLYPMAT, this function applies the projection on sources 
%     - The flow is computed via getVelocityField using the Stokeslet representation.
%     - Image systems are not included in this implementation.
%
%   See also: applyPmat, getCoarseSource, get_long_range_precond,
%   getVelocityField
%
% Anna Broms, Oct 19, 2025



%compute velocities using these source strengths

%remember that everything so far is implemented without images in mind
proj_vel = getVelocityField(rvec_in,rvec_out,vel(1:end/2),vel(end/2+1:end),[],[],[],[],[],[], []);

lambda = getCoarseSource(proj_vel,Rinv,Nx,Ny,Mx,Zi,Yi,opt.db,opt.P,opt.N_c,opt.a_c);
%lambda = An*Rinv*(Am'*proj_vel);

lambda_fine = vel-lambda;


end