function [lam_c,lam_self,lam_f,lam_e,u_corr,pair_qv_nonp, ...
    lam_c_nonp,lam_self_nonp,lam_f_nonp,lam_e_nonp] = ...
    transform_lap_peanut_big_sparse(tau,geom,basis)
%TRANSFORM_LAP_PEANUT_BIG_SPARSE Apply prebuilt sparse pair maps.

if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('transform_lap_peanut_big_sparse:MissingBigSparse', ...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('transform_lap_peanut_big_sparse:UnsupportedGrid', ...
        'The Laplace big-sparse transform is only supported on the solve grid.');
end

opt = geom.opt;
q = geom.q(:);
rvec_out = geom.rvec_out;
U = basis.U;
Y = basis.Y;

P = numel(q);
N_c = opt.N_c;
N_large = numel(rvec_out)/P;
n_coarse = P*N_c;
project_charge = logical(getOptField(opt,'project_charge',false));

lam_c = zeros(n_coarse,1);
lam_c_nonp = zeros(n_coarse,1);

for i = 1:P
    block = (i-1)*N_large+1:i*N_large;
    tau_i = tau(block);
    lam_i_nonp = Y{1}*(U{1}*tau_i);
    lam_i = projectChargeMode(lam_i_nonp,project_charge);
    idx = (i-1)*N_c+1:i*N_c;
    lam_c(idx) = lam_i;
    lam_c_nonp(idx) = lam_i_nonp;
end

lam_self = lam_c;
lam_self_nonp = lam_c_nonp;

[pair_nonp,pair_qv] = applyPairMaps(basis.big_sparse,lam_self);
pair_proj = applyDensePairProjector(pair_nonp,basis.big_sparse.P_pair);
corr = scatterSourceCorrections(basis.big_sparse.source_scatter_rows, ...
    pair_proj,pair_nonp,2*N_c,n_coarse);

lam_c = lam_self + corr(1:n_coarse);
lam_c_nonp = lam_self_nonp + corr(n_coarse+1:2*n_coarse);

pair_qv_nonp = scatterPairQv(basis.big_sparse.qv_scatter_rows,pair_qv,P);
u_corr = basis.big_sparse.M_u_cross*lam_self - ...
    basis.big_sparse.M_u_peanut*pair_proj;
if isfield(basis.big_sparse,'M_u_qv') && ~isempty(basis.big_sparse.M_u_qv)
    u_corr = u_corr - basis.big_sparse.M_u_qv*pair_qv;
end

lam_f = [];
lam_e = [];
lam_f_nonp = [];
lam_e_nonp = [];
end

function [pair_nonp,pair_qv] = applyPairMaps(big_sparse,lam_self)
if isfield(big_sparse,'pair_map_mode') && ...
        strcmp(big_sparse.pair_map_mode,'rotated_factored')
    pair_rhs = rotatePairInputs(lam_self,big_sparse);
    pair_nonp_canon = big_sparse.M_pair_nonp_canon*pair_rhs;
    pair_nonp = rotatePairOutputs(pair_nonp_canon,big_sparse);
    pair_qv = big_sparse.M_pair_qv_canon*pair_rhs;
else
    pair_nonp = big_sparse.M_pair_nonp*lam_self;
    pair_qv = big_sparse.M_pair_qv*lam_self;
end
end

function pair_rhs = rotatePairInputs(lam_self,big_sparse)
N_c = big_sparse.N_c;
pairs = big_sparse.pairs;
meta = big_sparse.pair_rotation_meta;
n_pairs = size(pairs,1);
pair_rhs = zeros(n_pairs*2*N_c,1);

for row = 1:n_pairs
    i = pairs(row,1);
    j = pairs(row,2);
    idx_i = (i-1)*N_c+1:i*N_c;
    idx_j = (j-1)*N_c+1:j*N_c;
    rhs = [lam_self(idx_i) lam_self(idx_j)];
    if ~isIdentityRotation(meta(row).rot)
        rhs = rotateUniformCircleData(rhs,[],meta(row).phase_c);
    end
    pair_idx = (row-1)*2*N_c+1:row*2*N_c;
    pair_rhs(pair_idx) = rhs(:);
end
end

function pair_nonp = rotatePairOutputs(pair_nonp_canon,big_sparse)
N_c = big_sparse.N_c;
meta = big_sparse.pair_rotation_meta;
n_pairs = numel(pair_nonp_canon)/(2*N_c);
pair_nonp = zeros(size(pair_nonp_canon));

for row = 1:n_pairs
    pair_idx = (row-1)*2*N_c+1:row*2*N_c;
    tau_pair = reshape(pair_nonp_canon(pair_idx),N_c,2);
    if ~isIdentityRotation(meta(row).rot)
        tau_pair = rotateUniformCircleData(tau_pair,[], ...
            meta(row).phase_c_inv);
    end
    pair_nonp(pair_idx) = tau_pair(:);
end
end

function tf = isIdentityRotation(rot)
tf = abs(rot - 1) <= 100*eps(max(1,abs(rot)));
end

function pair_proj = applyDensePairProjector(pair_nonp,P_pair)
if isempty(pair_nonp)
    pair_proj = pair_nonp;
    return
end
pair_rows = size(P_pair,1);
n_pairs = numel(pair_nonp)/pair_rows;
pair_proj = P_pair*reshape(pair_nonp,pair_rows,n_pairs);
pair_proj = pair_proj(:);
end

function corr = scatterSourceCorrections(rows,pair_proj,pair_nonp, ...
    pair_rows,n_coarse)
if isempty(rows)
    corr = zeros(2*n_coarse,1);
    return
end
n_pairs = numel(pair_nonp)/pair_rows;
pair_proj = reshape(pair_proj,pair_rows,n_pairs);
pair_nonp = reshape(pair_nonp,pair_rows,n_pairs);
vals = reshape([pair_proj; pair_nonp],[],1);
corr = accumarray(rows,vals,[2*n_coarse,1]);
end

function pair_qv_nonp = scatterPairQv(rows,pair_qv,P)
if isempty(rows)
    pair_qv_nonp = zeros(P,1);
else
    pair_qv_nonp = accumarray(rows,pair_qv(:),[P,1]);
end
end

function lam_out = projectChargeMode(lam_in,project_charge)
if ~project_charge || isempty(lam_in)
    lam_out = lam_in;
    return
end
n = numel(lam_in);
lam_out = lam_in - (sum(lam_in)/n);
end
