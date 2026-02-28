function lambda_fine = applyQmat(vel,rvec_in,rvec_out,Sinv,q,Ny,Mx,Zi,Yi,opt)
%APPLYQMAT  Apply long-range preconditioning projection matrix on sources
%
%   lambda_fine = APPLYQMAT(vel, rvec_in, rvec_out, Sinv, Zi, Yi, opt)
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
%     Sinv      - (Pk×Pk) inverse of the coarse interaction matrix S,
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

%mask out elements here
if opt.mask
    Nc = opt.N_c;
    a = opt.a_c;
    P = opt.P;
    vel_k = zeros(P*Nc*a*2,1);
    proj_vel = zeros(P*Nc*a*2,1);
    for k = 1:P
        vel_k(:) = 0;
        lambda_x = vel(Nc*(k-1)+1:Nc*k);
        lambda_y = vel(Nc*(k-1)+P*Nc+1:Nc*k+P*Nc);
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
    proj_vel = getVelocityField(rvec_in,rvec_out,vel(1:end/2),vel(end/2+1:end),[],[],[],[],[],[], []);
end



lambda = getCoarseSource(proj_vel,Sinv,Zi,Yi,opt.db,opt.P,opt.N_c,opt.a_c);

lambda_fine = vel-lambda;


end