function [Ucross,Ecolloc] = buildStokesCoarsePairDense( ...
    q,rbase_in_c,rout_base_c,pairs,row)
%BUILDSTOKESCOARSEPAIRDENSE Dense maps for one actual close-pair geometry.

i = pairs(row,1);
j = pairs(row,2);

rin_pair_c = [q(i)+rbase_in_c(:); q(j)+rbase_in_c(:)];
rout_pair = [q(i)+rout_base_c(:); q(j)+rout_base_c(:)];

Ecolloc = stokSLPmat(rin_pair_c,rout_pair,1);
Ucross = buildStokesCrossPairVelocityMap(Ecolloc,numel(rbase_in_c), ...
    numel(rout_base_c));
end
