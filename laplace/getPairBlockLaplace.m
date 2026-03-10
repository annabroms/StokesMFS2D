function [U,Y,A] = getPairBlockLaplace(rin_pair,rout_f)
%GETPAIRBLOCKLAPLACE Build pseudoinverse factors for scalar pair block.
%
% Syntax:
%   [U,Y,A] = getPairBlockLaplace(rin_pair,rout_f)
%
% Anna Broms, Mar 2026

A = lapSLPmat(rin_pair,rout_f);
[Y,U] = getPseudoFactors(A,1e-14,0);

end
