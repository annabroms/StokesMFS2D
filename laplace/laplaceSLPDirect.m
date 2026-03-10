function u = laplaceSLPDirect(rsrc,rtar,sigma)
%LAPLACESLPDIRECT Evaluate 2D Laplace single-layer potential by direct summation.
%
% Syntax:
%   u = laplaceSLPDirect(rsrc,rtar,sigma)
%
% Inputs:
%   rsrc  - Complex source locations.
%   rtar  - Complex target locations.
%   sigma - Source strengths, size numel(rsrc) x nd.
%
% Output:
%   u     - Potential values at targets, size numel(rtar) x nd.
%
% Notes:
%   Uses kernel -(1/(2*pi))*log(|x-y|). 
%
% Anna Broms, Mar 2026

rsrc = rsrc(:);
rtar = rtar(:);

ns = numel(rsrc);
nt = numel(rtar);

assert(size(sigma,1)==ns,'sigma must have one row per source location.');

nd = size(sigma,2);
u = zeros(nt,nd);

if ns==0 || nt==0
    return;
end


for k = 1:ns
    rk = abs(rtar-rsrc(k));   
    u = u - (1/(2*pi))*log(rk).*sigma(k,:);    
end

end
