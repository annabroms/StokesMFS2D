function lambda_fine = applyQmat_mob(vel,rvec_in,rvec_out,L,Lr,Sinv,Zi,Yi,opt)
%APPLYQMAT_MOB  Apply long-range preconditioning projection matrix on
%sources for the mobility setting
%
%   lambda_fine = APPLYQMAT_mob(vel, rvec_in, rvec_out, Sinv, Zi, Yi, opt)
%
%   Applies the projection operator 
%%   
%       Q = I - (Z * (S \ I) * Y') * G
%
%   to the solution vector lambda, with Z = blkdiag(Zi) and Y = blkdiag(Yi)
%
%   INPUTS:
%     lambda    - 3PN×1 vector of proxy source strengths 
%     rvec_in   - PN complex valued array of source points for evaluating the proxy flow.
%     rvec_out  - PM complex valued array of target surface points where the flow is evaluated.
%     Sinv      - (Pk×Pk) inverse of the coarse interaction matrix R,
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
%   See also: applyPmat_mob, getCoarseSource, get_long_range_precond_mob,
%   getVelocityField
%
% Anna Broms, Oct 22, 2025



%compute velocities using these source strengths
P = opt.P;
Nc = opt.N_c; 
a = opt.a_c;

proj_L = zeros(P*Nc*2,1);
Lr_vel = zeros(P*Nc*a*2,1);
for k = 1:opt.P
    vel_k = [vel(Nc*(k-1)+1:Nc*k); vel(P*Nc+Nc*(k-1)+1:Nc*k+P*Nc)];
    proj_vel_k = vel_k-L*vel_k; %project by I-L
    proj_L(Nc*(k-1)+1:Nc*k) = proj_vel_k(1:end/2);
    proj_L(Nc*(k-1)+P*Nc+1:Nc*k+P*Nc) = proj_vel_k(end/2+1:end);

    Lr_vel_k = Lr*vel_k; % Lr lambda
    Lr_vel(Nc*a*(k-1)+1:Nc*a*k) = Lr_vel_k(1:end/2);
    Lr_vel(Nc*a*(k-1)+P*Nc*a+1:Nc*a*k+P*Nc*a) = Lr_vel_k(end/2+1:end);


end

%remember that everything so far is implemented without images in mind
proj_vel = getVelocityField(rvec_in,rvec_out,proj_L(1:end/2),proj_L(end/2+1:end));

tot_vel = proj_vel+Lr_vel;

lambda = getCoarseSource(tot_vel,Sinv,Zi,Yi,opt.db,opt.P,opt.N_c,opt.a_c);

lambda_fine = vel-lambda;


end