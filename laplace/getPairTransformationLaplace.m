function [rvec_in,coarse_ind,lam_all,u_corr,lam_c,lam_f,lam_e] = getPairTransformationLaplace(tau,geom,basis)
%GETPAIRTRANSFORMATIONLAPLACE Map coarse boundary data to coarse/fine Laplace sources.
%
% Syntax:
%   [rvec_in,coarse_ind,lam_all,u_corr,lam_c,lam_f,lam_e] = ...
%       getPairTransformationLaplace(tau,geom,basis)
%
% Anna Broms, Mar 2026

rvec_out = geom.rvec_out;
q = geom.q;
pairs = geom.pairs;

rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
rimage_vec = geom.rimage_vec;
refine = geom.refine;
opt = geom.opt;

U = basis.U;
Y = basis.Y;
Upf = basis.Upf;
Ypf = basis.Ypf;

P = length(q);
N_c = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P;
precomp = opt.precomp;

lam_c = zeros(P*N_c,1);
rvec_in_coarse = zeros(P*N_c,1);
coarse_ind = cell(P,1);

lam_f = cell(P,1);
lam_e = cell(P,1);
lam_e_chunks = repmat({cell(0,1)},P,1);
rimage_k_chunks = repmat({cell(0,1)},P,1);
for k = 1:P
    lam_f{k} = zeros(N_f,1);
    lam_e{k} = zeros(0,1);
end

u_corr = zeros(P*N_large,1);

% Phase 1: one-body mapping from boundary data to coarse source strengths.
for i = 1:P
    tau_i = tau((i-1)*N_large+1:i*N_large);
    step1 = U{1}*tau_i;
    lam_i = Y{1}*step1;

    idx = (i-1)*N_c+1:i*N_c;
    coarse_ind{i} = idx;
    lam_c(idx) = lam_i;
    rvec_in_coarse(idx) = q(i)+rbase_in_c;
end

% Phase 2: pair corrections.
for row = 1:size(pairs,1)
    i = pairs(row,1);
    p2 = pairs(row,2);

    lam_i = lam_c(coarse_ind{i});
    lam_p2 = lam_c(coarse_ind{p2});

    if precomp
        rhs_pair = [lam_i; lam_p2];
    else
        rout_fine_other2 = getFineOther(opt.a_f,opt.N_f,refine,q,i,p2,opt.alpha);
        R2 = -lapSLPmat(rbase_in_c+q(i),rout_fine_other2)*lam_i;

        rout_fine_other1 = getFineOther(opt.a_f,opt.N_f,refine,q,p2,i,opt.alpha);
        R1 = -lapSLPmat(rbase_in_c+q(p2),rout_fine_other1)*lam_p2;

        rhs_pair = [R1; R2];
    end

    pair_mapped = Upf{i,p2}*rhs_pair;
    beta_tot = Ypf{i,p2}*pair_mapped;

    rimage_i = rimage_vec{i,p2};
    rimage_p2 = rimage_vec{p2,i};
    im_i = length(rimage_i);
    im_p2 = length(rimage_p2);

    s_i = 1:N_f;
    e_i = N_f+1:N_f+im_i;
    s_p2 = N_f+im_i+1:2*N_f+im_i;
    e_p2 = 2*N_f+im_i+1:2*N_f+im_i+im_p2;

    lam_f{i} = lam_f{i} + beta_tot(s_i);
    lam_f{p2} = lam_f{p2} + beta_tot(s_p2);
    lam_e_chunks{i}{end+1,1} = beta_tot(e_i);
    lam_e_chunks{p2}{end+1,1} = beta_tot(e_p2);
    rimage_k_chunks{i}{end+1,1} = rimage_i;
    rimage_k_chunks{p2}{end+1,1} = rimage_p2;

    % Pair-local flow contribution.
    block_i = (i-1)*N_large+1:i*N_large;
    block_p2 = (p2-1)*N_large+1:p2*N_large;
    rout_pair = [rvec_out(block_i); rvec_out(block_p2)];
    rin_pair = [q(i)+rbase_in_f; rimage_i; q(p2)+rbase_in_f; rimage_p2];

    u_pair = laplaceSingleLayerField(rin_pair,rout_pair,beta_tot,false);

    pair_ind = [block_i block_p2]';
    u_corr(pair_ind) = u_corr(pair_ind) + u_pair;

    % Add coarse cross interactions, so pair blocks can be identity-corrected in matvec.
    u2 = lapSLPmat(rbase_in_c+q(i),rvec_out(block_p2))*lam_i;
    u1 = lapSLPmat(rbase_in_c+q(p2),rvec_out(block_i))*lam_p2;

    u_corr(block_i) = u_corr(block_i) + u1;
    u_corr(block_p2) = u_corr(block_p2) + u2;
end

% Assemble all source points and strengths.
has_neigh = sort(unique(pairs(:)));

r_chunks = {rvec_in_coarse};
l_chunks = {lam_c};

for ii = 1:numel(has_neigh)
    k = has_neigh(ii);
    if isempty(lam_e_chunks{k})
        lam_e{k} = zeros(0,1);
        rimage_k = zeros(0,1);
    else
        lam_e{k} = vertcat(lam_e_chunks{k}{:});
        rimage_k = vertcat(rimage_k_chunks{k}{:});
    end

    r_chunks{end+1,1} = q(k)+rbase_in_f; 
    l_chunks{end+1,1} = lam_f{k}; 
    r_chunks{end+1,1} = rimage_k; 
    l_chunks{end+1,1} = lam_e{k}; 
end

rvec_in = vertcat(r_chunks{:});
lam_all = vertcat(l_chunks{:});

end
