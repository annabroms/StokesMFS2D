function group = buildLaplacePairGroup(group_id,row,sep,q,rbase_in_c, ...
    rbase_in_f,rimage_vec,refine,pairs,opt,rout_base_f,project_charge, ...
    use_canonical)
%BUILDLAPLACEPAIRGROUP Build one Laplace pair/group basis payload.

if nargin < 13
    use_canonical = true;
end

q = q(:);
i = pairs(row,1);
j = pairs(row,2);
mid = 0.5*(q(i)+q(j));
delta = q(j)-q(i);
if abs(delta) == 0
    rot = 1;
else
    rot = delta/abs(delta);
end
R = opt.rad;

if use_canonical
    q_pair = [-sep/2; sep/2];
    rimage_i = map_points_to_canonical(rimage_vec{i,j},mid,rot);
    rimage_j = map_points_to_canonical(rimage_vec{j,i},mid,rot);
    refine_i = map_points_to_canonical(refine{i,j},mid,rot);
    refine_j = map_points_to_canonical(refine{j,i},mid,rot);
else
    q_pair = [q(i); q(j)];
    rimage_i = rimage_vec{i,j};
    rimage_j = rimage_vec{j,i};
    refine_i = refine{i,j};
    refine_j = refine{j,i};
end

rout_f = [q_pair(1)+rout_base_f; refine_i; q_pair(2)+rout_base_f; refine_j];
rin_pair_f = [q_pair(1)+rbase_in_f; rimage_i; ...
              q_pair(2)+rbase_in_f; rimage_j];
nsrc_f_i = numel(rbase_in_f) + numel(rimage_i);
nsrc_f_j = numel(rbase_in_f) + numel(rimage_j);
ntar_i = numel(rout_base_f) + numel(refine_i);
ntar_j = numel(rout_base_f) + numel(refine_j);

proj_pair = struct('project_charge',project_charge, ...
    'nsrc',[nsrc_f_i nsrc_f_j], ...
    'ntar',[ntar_i ntar_j]);
[Uf_pair,Yf_pair] = getPairBlockLaplace(rin_pair_f,rout_f,proj_pair);

Npair = evaluateCoarseOnPairLaplace(q_pair,rbase_in_c,rout_f);
Upf = -Uf_pair'*Npair;
Ypf = Yf_pair;

DC = [];
YC = [];
C = [];
QV = [];

if opt.N_peanut
    rout_peanut = createPeanut(q_pair(1),q_pair(2),opt.N_peanut,0,R);
    rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
    proj_peanut = struct('project_charge',project_charge, ...
        'nsrc_c',[numel(rbase_in_c) numel(rbase_in_c)], ...
        'nsrc_f',[nsrc_f_i nsrc_f_j]);
    [DC,YC] = getPeanutBlockLaplace(rin_pair_c,rin_pair_f, ...
        rout_peanut,proj_peanut);

    if opt.cmap
        C = -YC*(DC*Yf_pair*(Uf_pair'*Npair));
        QV = getPairChargeSumMap(nsrc_f_i,nsrc_f_j)*Ypf*Upf;
        if opt.compress_cmap
            [U,S,V] = svd(C);
            S = diag(S);
            ra = sum(S>max(S)*opt.cmap_tol);
            fprintf('Rank of coarse-coarse map is chosen to %u\n',ra);
            S = S(1:ra);
            C = U(:,1:ra)*diag(S)*V(:,1:ra)';
        end
    end
end

group = struct('group_id',group_id,'sep',sep,'q_pair',q_pair, ...
    'rimage_canon',{{rimage_i; rimage_j}}, ...
    'refine_canon',{{refine_i; refine_j}}, ...
    'Upf',Upf,'Ypf',Ypf,'DC',DC,'YC',YC,'Cmap',C, ...
    'Cmap_QV',QV,'nout_f',numel(rout_base_f), ...
    'nsrc_f',[nsrc_f_i nsrc_f_j],'ntar_f',[ntar_i ntar_j], ...
    'rep_pair',[i j]);
end

function Kq_pair = getPairChargeSumMap(n1,n2)
Kq_pair = [ones(1,n1) zeros(1,n2); zeros(1,n1) ones(1,n2)];
end

function z_canon = map_points_to_canonical(z,mid,rot)
z = z(:);
if isempty(z)
    z_canon = zeros(0,1);
    return
end

z_canon = conj(rot)*(z-mid);
end
