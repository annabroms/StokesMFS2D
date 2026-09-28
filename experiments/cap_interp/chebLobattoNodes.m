function s = chebLobattoNodes(smin,smax,q)
%CHEBLOBATTONODES Chebyshev-Lobatto (Clenshaw-Curtis) nodes on [smin,smax].
%
% s = chebLobattoNodes(smin,smax,q) returns the q nodes
%   s_j = (smin+smax)/2 + (smax-smin)/2*cos(j*pi/(q-1)),  j=0,...,q-1
% in that order (s(1)=smax, s(end)=smin). This order matches the standard
% barycentric weight formula used by chebBarycentricWeights/
% evalMatrixChebBary; do not resort the nodes.
%
% For q of the form 2^k+1 (e.g. 9,17,33,65), these grids are nested: every
% node of a level q also appears (to machine precision) in the next level
% 2*(q-1)+1.
%
% Anna Broms, Sep 2026

j = (0:q-1)';
s = (smin+smax)/2 + (smax-smin)/2*cos(j*pi/(q-1));

end
