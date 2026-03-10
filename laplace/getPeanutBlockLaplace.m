function [DC,Y] = getPeanutBlockLaplace(rin_pair_c,rin_pair_f,rout_peanut)
%GETPEANUTBLOCKLAPLACE Factorization for scalar peanut compression.
%
% Syntax:
%   [DC,Y] = getPeanutBlockLaplace(rin_pair_c,rin_pair_f,rout_peanut)

Nf = lapSLPmat(rin_pair_f,rout_peanut);
Nc = lapSLPmat(rin_pair_c,rout_peanut);

[Y,U] = getPseudoFactors(Nc,1e-14,0);
DC = U'*Nf;

end
