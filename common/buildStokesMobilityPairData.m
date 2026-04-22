function [pair,debug_data] = buildStokesMobilityPairData( ...
    q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,Lc,row, ...
    debug,svd_opts,use_canonical,payload_mode)
%BUILDSTOKESMOBILITYPAIRDATA Build one Stokes mobility pair payload.
%
% This helper contains the pair solve used by getPairBasisStokes and by the
% streamed Stokes big-sparse builder. The caller decides whether to keep the
% full fine-source payload or only the compressed maps.

if nargin < 10 || isempty(debug)
    debug = false;
end
if nargin < 11 || isempty(svd_opts)
    svd_opts = struct();
end
if nargin < 12 || isempty(use_canonical)
    use_canonical = false;
end
if nargin < 13 || isempty(payload_mode)
    payload_mode = 'full';
end
store_full_pair_payload = strcmp(payload_mode,'full');

q = q(:);
rbase_in_c = rbase_in_c(:);
rbase_in_f = rbase_in_f(:);

i = pairs(row,1);
j = pairs(row,2);
pair = init_pair_payload(i,j);

delta = q(j)-q(i);
sep = abs(delta);
if sep == 0
    rot = 1;
else
    rot = delta/sep;
end
mid = 0.5*(q(i)+q(j));

if use_canonical
    q_pair = [-sep/2; sep/2];
    rimage_i = map_points_to_canonical(rimage_pairs{i,j},mid,rot);
    rimage_j = map_points_to_canonical(rimage_pairs{j,i},mid,rot);
    refine_i = map_points_to_canonical(refine{i,j},mid,rot);
    refine_j = map_points_to_canonical(refine{j,i},mid,rot);
else
    q_pair = [q(i); q(j)];
    rimage_i = rimage_pairs{i,j};
    rimage_j = rimage_pairs{j,i};
    refine_i = refine{i,j};
    refine_j = refine{j,i};
end

pair.q_pair = q_pair;
pair.rimage_canon = {rimage_i; rimage_j};
pair.refine_canon = {refine_i; refine_j};

project_force = logical(getOptField(opt,'project_force',false));

N_c = opt.N_c;
a_c = getOptField(opt,'a_c',1.2);
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';
rout_base_c = cos(tout_c)+1i*sin(tout_c);

nout = ceil(opt.a_f*opt.N_f);
t = linspace(0,2*pi,nout+1);
t = t(1:end-1)';
rout_base = cos(t)+1i*sin(t);

rin_1_f = q_pair(1)+rbase_in_f;
rin_2_f = q_pair(2)+rbase_in_f;
rout_f = [q_pair(1)+rout_base; refine_i; q_pair(2)+rout_base; refine_j];
rin_pair = [rin_1_f; rimage_i; rin_2_f; rimage_j];
if store_full_pair_payload
    pair.rin_pair = rin_pair;
end

svd_pair = svd_opts;
if logical(getOptField(svd_opts,'left_weight',false))
    row_weights_1 = getPeriodicCurveWeights( ...
        [q_pair(1)+rout_base; refine_i],q_pair(1));
    row_weights_2 = getPeriodicCurveWeights( ...
        [q_pair(2)+rout_base; refine_j],q_pair(2));
    svd_pair.row_weights = [row_weights_1; row_weights_2];
end

if debug
    plot_pair_geometry(i,j,q_pair,rin_1_f,rin_2_f,rout_base, ...
        refine_i,refine_j,rimage_i,rimage_j);
end

if ~isempty(Lc)
    Lc_pair = getILpair(Lc);
else
    Lc_pair = [];
end

if project_force || logical(getOptField(opt,'cmap',false))
    Kf1 = getKmat2D(rin_pair(1:end/2),q_pair(1));
    Kf2 = getKmat2D(rin_pair(end/2+1:end),q_pair(2));
else
    Kf1 = [];
    Kf2 = [];
end

if project_force
    B1 = getKmat2D([q_pair(1)+rout_base; refine_i],q_pair(1));
    B2 = getKmat2D([q_pair(2)+rout_base; refine_j],q_pair(2));
    Lr_pair = getLrPair(B1,B2,Kf1,Kf2);
    Lf_pair = getLfPair(Kf1,Kf2);
else
    Lr_pair = [];
    Lf_pair = [];
end

[Uf_pair,Yf_pair] = getPairBlockStokes(rin_pair,rout_f,Lf_pair, ...
    Lr_pair,svd_pair);
Npair = evaluateCoarseOnPair(q_pair,rbase_in_c,rout_f);
Upf = -Uf_pair'*Npair;

DC = [];
YC = [];
if getOptField(opt,'N_peanut',0)
    rout_peanut = createPeanut(q_pair(1),q_pair(2),opt.N_peanut,0);
    rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
    [DC,YC] = getPeanutBlockStokes(rin_pair_c,rin_pair,rout_peanut, ...
        Lc_pair,Lf_pair,svd_opts);
    if logical(getOptField(opt,'cmap',false))
        pair.Cmap = -YC*(DC*Yf_pair*(Uf_pair'*Npair));
    end
else
    rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
end

if logical(getOptField(opt,'cmap',false))
    Kft_pair = getKftPair(Kf1,Kf2);
    pair.Cmap_FU = -Kft_pair*Yf_pair*(Uf_pair'*Npair);
end

if store_full_pair_payload
    pair.Upf = Upf;
    pair.Ypf = Yf_pair;
    pair.Lf_pair = Lf_pair;
    pair.Lc_pair = Lc_pair;
    pair.DC = DC;
    pair.YC = YC;
end

if store_full_pair_payload && project_force
    rout_pair_c = [q_pair(1)+rout_base_c; q_pair(2)+rout_base_c];
    Ppair = eye(size(Lf_pair)) - Lf_pair;
    pair.Upair_colloc = stokSLPmat(rin_pair,rout_pair_c,1) * ...
        Ppair * Yf_pair * Upf;
    pair.Ecolloc = stokSLPmat(rin_pair_c,rout_pair_c,1);
    pair.Ucross_colloc = build_cross_pair_velocity_map( ...
        pair.Ecolloc,N_c,numel(rout_base_c));
end

debug_data = struct();
debug_data.q_pair = q_pair;
debug_data.rin_pair = rin_pair;
debug_data.Upf = Upf;
debug_data.Yf_pair = Yf_pair;
debug_data.Lf_pair = Lf_pair;
debug_data.Kf1 = Kf1;
debug_data.Kf2 = Kf2;
debug_data.Lc_pair = Lc_pair;
debug_data.DC = DC;
debug_data.YC = YC;
debug_data.rout_f = rout_f;
debug_data.Npair = Npair;
end

function pair = init_pair_payload(i,j)
pair = struct();
pair.rep_pair = [i j];
pair.q_pair = [];
pair.rimage_canon = {[]; []};
pair.refine_canon = {[]; []};
pair.rin_pair = [];
pair.Upf = [];
pair.Ypf = [];
pair.Lf_pair = [];
pair.Lc_pair = [];
pair.DC = [];
pair.YC = [];
pair.Cmap = [];
pair.Cmap_FU = [];
pair.Upair_colloc = [];
pair.Ucross_colloc = [];
pair.Ecolloc = [];
end

function plot_pair_geometry(i,j,q_pair,rin_1_f,rin_2_f,rout_base, ...
    refine_i,refine_j,rimage_i,rimage_j)
figure(801);
clf;
plot(real(rin_1_f),imag(rin_1_f),'r.','MarkerSize',10);
hold on;
plot(real(rin_2_f),imag(rin_2_f),'b.','MarkerSize',10);
plot(real(q_pair(1)+rout_base),imag(q_pair(1)+rout_base),'ro', ...
    'MarkerSize',4);
plot(real(q_pair(2)+rout_base),imag(q_pair(2)+rout_base),'bo', ...
    'MarkerSize',4);
plot(real(refine_i),imag(refine_i),'r+','MarkerSize',6);
plot(real(refine_j),imag(refine_j),'b+','MarkerSize',6);
if ~isempty(rimage_i)
    plot(real(rimage_i),imag(rimage_i),'ks','MarkerSize',5);
end
if ~isempty(rimage_j)
    plot(real(rimage_j),imag(rimage_j),'kd','MarkerSize',5);
end
axis equal;
grid on;
title(sprintf('getPairBasisStokes pair (%d,%d)',i,j), ...
    'Interpreter','none');
drawnow;
end

function Ucross = build_cross_pair_velocity_map(Epair,N_src,N_tgt)
Ucross = zeros(size(Epair));

tgt_i_x = 1:N_tgt;
tgt_j_x = N_tgt+1:2*N_tgt;
tgt_i_y = 2*N_tgt+1:3*N_tgt;
tgt_j_y = 3*N_tgt+1:4*N_tgt;

src_i_x = 1:N_src;
src_j_x = N_src+1:2*N_src;
src_i_y = 2*N_src+1:3*N_src;
src_j_y = 3*N_src+1:4*N_src;

Ucross(tgt_i_x,[src_j_x src_j_y]) = Epair(tgt_i_x,[src_j_x src_j_y]);
Ucross(tgt_j_x,[src_i_x src_i_y]) = Epair(tgt_j_x,[src_i_x src_i_y]);
Ucross(tgt_i_y,[src_j_x src_j_y]) = Epair(tgt_i_y,[src_j_x src_j_y]);
Ucross(tgt_j_y,[src_i_x src_i_y]) = Epair(tgt_j_y,[src_i_x src_i_y]);

Ucross = -Ucross;
end

function z_canon = map_points_to_canonical(z,mid,rot)
z = z(:);
if isempty(z)
    z_canon = zeros(0,1);
    return
end
z_canon = conj(rot)*(z-mid);
end
