function w = chebBarycentricWeights(q)
%CHEBBARYCENTRICWEIGHTS 2nd-kind barycentric weights for chebLobattoNodes(...,q).
%
% w = chebBarycentricWeights(q) returns w_j = (-1)^j, j=0,...,q-1, with the
% two endpoint weights halved. Must be paired with nodes produced by
% chebLobattoNodes (same order, not resorted).
%
% Anna Broms, Sep 2026

j = (0:q-1)';
w = (-1).^j;
w([1 end]) = w([1 end])/2;

end
