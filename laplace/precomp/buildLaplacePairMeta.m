function meta = buildLaplacePairMeta(q,pairs,nc,nf,opt)
%BUILDLAPLACEPAIRMETA Build per-actual-pair rotation and grouping metadata.

q = q(:);
total_pairs = size(pairs,1);
meta = repmat(struct('i',[],'j',[],'group_id',[],'sep',[], ...
    'mid',[],'rot',[],'phase_c',[],'phase_c_inv',[], ...
    'phase_f',[],'phase_f_inv',[], ...
    'Ucross_colloc_actual',[],'Ec_colloc_actual',[], ...
    'Lr_colloc_actual',[]),total_pairs,1);

for row = 1:total_pairs
    i = pairs(row,1);
    j = pairs(row,2);
    delta = q(j)-q(i);
    sep = abs(delta);
    if sep == 0
        rot = 1;
    else
        rot = delta/sep;
    end

    meta(row).i = i;
    meta(row).j = j;
    meta(row).sep = sep;
    meta(row).mid = 0.5*(q(i)+q(j));
    meta(row).rot = rot;
    meta(row).phase_c = getUniformCircleRotationSpec(nc,rot,opt);
    meta(row).phase_c_inv = invertUniformCircleRotationSpec(meta(row).phase_c);
    meta(row).phase_f = getUniformCircleRotationSpec(nf,rot,opt);
    meta(row).phase_f_inv = invertUniformCircleRotationSpec(meta(row).phase_f);
end
end
