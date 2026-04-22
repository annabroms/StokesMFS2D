function [rvec_in,coarse_ind,lam_all,lam_c,lam_f,lam_e,lam_all_nonp,lam_c_nonp,lam_f_nonp,lam_e_nonp] = getPairTransformationLaplace(tau,geom,basis)
%GETPAIRTRANSFORMATIONLAPLACE Map coarse boundary data to interior coarse/fine Laplace sources.
%
% Syntax:
%   [rvec_in,coarse_ind,lam_all,lam_c,lam_f,lam_e,...
%    lam_all_nonp,lam_c_nonp,lam_f_nonp,lam_e_nonp] = ...
%       getPairTransformationLaplace(tau,geom,basis)
%
% Outputs:
%   rvec_in    - Complex valued array of all source points, sorted per
%                object
%   lam_*      - Projected source strengths (charge modes removed per body).
%   lam_*_nonp - Non-projected source strengths (needed for computation of voltages).
%
% See also: getPairBasisLaplace, matvec_lap_2B_enhanced, evaluateCoarseOnPairLaplace.
%
% Anna Broms, Mar 2026

rvec_out = geom.rvec_out;
q = geom.q;
pairs = geom.pairs;

rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
rimage_vec = geom.rimage_vec;
opt = geom.opt;

U = basis.U;
Y = basis.Y;
Upf = basis.Upf;
Ypf = basis.Ypf;
pair_cache = basis.pair_cache;
use_pair_cache = pair_cache.enabled;

if isfield(opt,'project_charge') && ~isempty(opt.project_charge)
    project_charge = logical(opt.project_charge);
else
    project_charge = false;
end

P = numel(q);
N_c = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P;

lam_c = zeros(P*N_c,1);
lam_c_nonp = zeros(P*N_c,1);
rvec_in_coarse = zeros(P*N_c,1);
coarse_ind = cell(P,1);

lam_f = cell(P,1);
lam_e = cell(P,1);
lam_f_nonp = cell(P,1);
lam_e_nonp = cell(P,1);

lam_e_chunks = repmat({cell(0,1)},P,1);
lam_e_nonp_chunks = repmat({cell(0,1)},P,1);
rimage_k_chunks = repmat({cell(0,1)},P,1);

for k = 1:P
    lam_f{k} = zeros(N_f,1);
    lam_e{k} = zeros(0,1);
    lam_f_nonp{k} = zeros(N_f,1);
    lam_e_nonp{k} = zeros(0,1);
end

% Phase 1: one-body mapping from boundary data to coarse source strengths.
for i = 1:P
    tau_i = tau((i-1)*N_large+1:i*N_large);
    step1 = U{1}*tau_i;
    lam_i_nonp = Y{1}*step1;

    if project_charge
        lam_i = lam_i_nonp - mean(lam_i_nonp);
    else
        lam_i = lam_i_nonp;
    end

    idx = (i-1)*N_c+1:i*N_c;
    coarse_ind{i} = idx;
    lam_c(idx) = lam_i;
    lam_c_nonp(idx) = lam_i_nonp;
    rvec_in_coarse(idx) = q(i)+rbase_in_c;
end

% Phase 2: pair corrections.
for row = 1:size(pairs,1)
    i = pairs(row,1);
    p2 = pairs(row,2);

    lam_i = lam_c(coarse_ind{i});
    lam_p2 = lam_c(coarse_ind{p2});

    if use_pair_cache
        pair = getLaplacePairInstance(pair_cache,row);
        phase_c = pair.meta.phase_c;
        phase_f_inv = pair.meta.phase_f_inv;
        group = pair.group;

        rhs_pair = rotateUniformCircleData([lam_i lam_p2],[],phase_c);
        %rhs_pair = [lam_i lam_p2]; %debug;
        rhs_pair = rhs_pair(:);
        

        pair_mapped = group.Upf*rhs_pair;
        beta_tot_nonp_local = group.Ypf*pair_mapped;

        rimage_i = pair.rimage_i;
        rimage_p2 = pair.rimage_j;
        im_i = numel(group.rimage_canon{1});
        im_p2 = numel(group.rimage_canon{2});
    else
        rhs_pair = [lam_i; lam_p2];
        pair_mapped = Upf{i,p2}*rhs_pair;
        beta_tot_nonp_local = Ypf{i,p2}*pair_mapped;

        rimage_i = rimage_vec{i,p2};
        rimage_p2 = rimage_vec{p2,i};
        im_i = length(rimage_i);
        im_p2 = length(rimage_p2);
    end

    s_i = 1:N_f;
    e_i = N_f+1:N_f+im_i;
    s_p2 = N_f+im_i+1:2*N_f+im_i;
    e_p2 = 2*N_f+im_i+1:2*N_f+im_i+im_p2;

    beta_i_nonp_local = [beta_tot_nonp_local(s_i); beta_tot_nonp_local(e_i)];
    beta_p2_nonp_local = [beta_tot_nonp_local(s_p2); beta_tot_nonp_local(e_p2)];

    if project_charge
        beta_i_local = beta_i_nonp_local - mean(beta_i_nonp_local);
        beta_p2_local = beta_p2_nonp_local - mean(beta_p2_nonp_local);
    else
        beta_i_local = beta_i_nonp_local;
        beta_p2_local = beta_p2_nonp_local;
    end

    beta_fine_pair = [beta_i_local(1:N_f) beta_p2_local(1:N_f) ...
        beta_i_nonp_local(1:N_f) beta_p2_nonp_local(1:N_f)];
    if use_pair_cache
        beta_fine_pair = rotateUniformCircleData(beta_fine_pair,[],phase_f_inv);
    end
    beta_i_fine = beta_fine_pair(:,1);
    beta_p2_fine = beta_fine_pair(:,2);
    beta_i_fine_nonp = beta_fine_pair(:,3);
    beta_p2_fine_nonp = beta_fine_pair(:,4);

    lam_f{i} = lam_f{i} + beta_i_fine;
    lam_f{p2} = lam_f{p2} + beta_p2_fine;
    lam_e_chunks{i}{end+1,1} = beta_i_local(N_f+1:end);
    lam_e_chunks{p2}{end+1,1} = beta_p2_local(N_f+1:end);

    lam_f_nonp{i} = lam_f_nonp{i} + beta_i_fine_nonp;
    lam_f_nonp{p2} = lam_f_nonp{p2} + beta_p2_fine_nonp;
    lam_e_nonp_chunks{i}{end+1,1} = beta_i_nonp_local(N_f+1:end);
    lam_e_nonp_chunks{p2}{end+1,1} = beta_p2_nonp_local(N_f+1:end);

    rimage_k_chunks{i}{end+1,1} = rimage_i;
    rimage_k_chunks{p2}{end+1,1} = rimage_p2;
end

% Assemble all source points and strengths.
has_neigh = sort(unique(pairs(:)));

r_chunks = {rvec_in_coarse};
l_chunks = {lam_c};
l_chunks_nonp = {lam_c_nonp};

for ii = 1:numel(has_neigh)
    k = has_neigh(ii);
    if isempty(lam_e_chunks{k})
        lam_e{k} = zeros(0,1);
        lam_e_nonp{k} = zeros(0,1);
        rimage_k = zeros(0,1);
    else
        lam_e{k} = vertcat(lam_e_chunks{k}{:});
        lam_e_nonp{k} = vertcat(lam_e_nonp_chunks{k}{:});
        rimage_k = vertcat(rimage_k_chunks{k}{:});
    end

    r_chunks{end+1,1} = q(k)+rbase_in_f;
    l_chunks{end+1,1} = lam_f{k};
    l_chunks_nonp{end+1,1} = lam_f_nonp{k};

    r_chunks{end+1,1} = rimage_k;
    l_chunks{end+1,1} = lam_e{k};
    l_chunks_nonp{end+1,1} = lam_e_nonp{k};
end

rvec_in = vertcat(r_chunks{:});
lam_all = vertcat(l_chunks{:});
lam_all_nonp = vertcat(l_chunks_nonp{:});

end
