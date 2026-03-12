function [udirect,vdirect] = stokSLPdirect(xsrc,ysrc,xtar,ytar,f1,f2,N)
%STOKSLPDIRECT  Evaluate 2D Stokeslet flow from discrete sources using direct summation.
%
%   [udirect, vdirect] = STOKSLPDIRECT(xsrc, ysrc, xtar, ytar, f1, f2, N)
%
%   Computes the velocity field (u,v) at a set of target points due to N
%   Stokeslet sources in 2D, located at (xsrc, ysrc) with force strengths (f1, f2).
%   The velocity is evaluated using the free-space 2D Stokeslet kernel:
%
%       u(x) = 1/(4π) ∑_k [ -½ log(r_k²) I + (r_k ⊗ r_k)/r_k² ] * f_k,
%
%   where r_k = x - x_k is the vector from source to target, and f_k is the force.
%
%   INPUTS:
%     xsrc    - N×1 vector of x-coordinates of sources.
%     ysrc    - N×1 vector of y-coordinates of sources.
%     xtar    - M×1 vector of x-coordinates of targets.
%     ytar    - M×1 vector of y-coordinates of targets.
%     f1      - N×nd array of x-components of source strengths.
%     f2      - N×nd array of y-components of source strengths.
%     N       - Number of source points (should match size of xsrc, ysrc, f1, f2).
%
%   OUTPUTS:
%     udirect - M×nd array of x-components of velocity at each target point.
%     vdirect - M×nd array of y-components of velocity at each target point.
%
%   NOTES:
%     - Sources and targets must not coincide (no self-interaction).
%     - The kernel is the fundamental solution to the 2D Stokes equations
%       in free space with viscosity μ = 1.
%     - Multiple force vectors (nd > 1) are supported (batched columns).
%
% Anna Broms, Oct 19, 2025

assert(numel(xsrc)==N,'Wrong number of x sources');
assert(numel(ysrc)==N,'Wrong number of y sources');
assert(size(f1,1)==N,'Wrong number of x source strengths');
assert(size(f2,1)==N,'Wrong number of y source strengths');

nd = size(f1,2);

udirect = zeros(length(ytar),nd);
vdirect = zeros(length(ytar),nd); 


for k = 1:N % loop over sources
     rx = xtar - xsrc(k);
     ry = ytar - ysrc(k);
     rho2 = rx.^2 + ry.^2;
     rdotf = rx.*f1(k,:) + ry.*f2(k,:);
     logrho2 = log(rho2);
     udirect = udirect-0.5*logrho2.*f1(k,:) + rdotf./rho2.*rx;
     vdirect = vdirect-0.5*logrho2.*f2(k,:) + rdotf./rho2.*ry;
end


    vdirect = vdirect/4/pi;
    udirect = udirect/4/pi; 
end
