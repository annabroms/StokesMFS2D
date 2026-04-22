function [Ucross,Ec,Lr] = buildLaplaceActualPairCollocFactors( ...
    meta,q,rbase_in_c,rout_base_c)
%BUILDLAPLACEACTUALPAIRCOLLOCFACTORS Dense solve-grid pair blocks.

q = q(:);
i = meta.i;
j = meta.j;
ntar_body = numel(rout_base_c);
nsrc_body = numel(rbase_in_c);

rout_pair = [q(i)+rout_base_c; q(j)+rout_base_c];
rin_pair_c = [q(i)+rbase_in_c; q(j)+rbase_in_c];

Ec = lapSLPmat(rin_pair_c,rout_pair);

Ucross = zeros(2*ntar_body,2*nsrc_body);
rows_i = 1:ntar_body;
rows_j = ntar_body+1:2*ntar_body;
cols_i = 1:nsrc_body;
cols_j = nsrc_body+1:2*nsrc_body;
Ucross(rows_i,cols_j) = -Ec(rows_i,cols_j);
Ucross(rows_j,cols_i) = -Ec(rows_j,cols_i);
Lr = [ones(ntar_body,1) zeros(ntar_body,1); ...
      zeros(ntar_body,1) ones(ntar_body,1)];
end
