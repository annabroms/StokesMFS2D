function [lam_c,lam_self,lam_f,lam_e,u_corr,lam_c_nonp,lam_self_nonp,lam_f_nonp,lam_e_nonp] = transform_lap_peanut(tau,geom,basis)
%TRANSFORM_LAP_PEANUT Map coarse boundary data to interior compressed/fine
%Laplace sources. Helper function to matvec_lap_peanut_enhanced
%
% Syntax:
%   [lam_c,lam_self,lam_f,lam_e,u_corr] = transform_lap_peanut(tau,geom,basis)
%   [lam_c,lam_self,lam_f,lam_e,u_corr,lam_c_nonp,lam_self_nonp, ...
%       lam_f_nonp,lam_e_nonp] = transform_lap_peanut(tau,geom,basis)
%
% Outputs:
%   lam_*      - Projected source strengths.
%   lam_*_nonp - Non-projected source strengths (used for Lr/voltage terms).
%
% See also: getPairBasisLaplace, getPeanutBlockLaplace, ...
%   matvec_lap_peanut_enhanced.
%
% Anna Broms, Mar 2026

rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
rimage_vec = geom.rimage_vec;
refine = geom.refine;
opt = geom.opt;
rvec_out = geom.rvec_out;
q = geom.q;
pairs = geom.pairs;
if isfield(geom,'rcheck') && ~isempty(geom.rcheck)
    rcheck_out = geom.rcheck;
else
    rcheck_out = rvec_out;
end

U = basis.U;
Y = basis.Y;
Upf = basis.Upf;
Ypf = basis.Ypf;
DC_all = basis.DC_all;
YC_all = basis.YC_all;
Cmap = basis.Cmap;
pair_cache = get_pair_cache(geom,basis);
use_pair_cache = ~isempty(pair_cache) && isfield(pair_cache,'enabled') && pair_cache.enabled;

P = numel(q);
N_c = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P;
N_check = length(rcheck_out)/P;
if isfield(opt,'project_charge') && ~isempty(opt.project_charge)
    project_charge = logical(opt.project_charge);
else
    project_charge = false;
end

precomp = opt.precomp;

lam_c = zeros(N_c*P,1);
lam_c_nonp = zeros(N_c*P,1);
lam_f = cell(P,1);
lam_e = cell(P,1);
lam_f_nonp = cell(P,1);
lam_e_nonp = cell(P,1);
lam_e_chunks = repmat({cell(0,1)},P,1);
lam_e_nonp_chunks = repmat({cell(0,1)},P,1);
for k = 1:P
    lam_f{k} = zeros(N_f,1);
    lam_e{k} = zeros(0,1);
    lam_f_nonp{k} = zeros(N_f,1);
    lam_e_nonp{k} = zeros(0,1);
end

u_corr = zeros(P*N_check,1);

% Phase 1: one-body coarse mapping.
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

% Phase 2: pair loop.
for row = 1:size(pairs,1)
    i = pairs(row,1);
    p2 = pairs(row,2);

    idx_i = (i-1)*N_c+1:i*N_c;
    idx_p2 = (p2-1)*N_c+1:p2*N_c;

    lam_i = lam_self(idx_i);
    lam_p2 = lam_self(idx_p2);

    if use_pair_cache
        pair = getLaplacePairInstance(pair_cache,row);
        rot = pair.meta.rot;
        group = pair.group;

        lam_i_local = rotateUniformCircleData(lam_i,rot);
        lam_p2_local = rotateUniformCircleData(lam_p2,rot);

        if precomp
            rhs_pair = [lam_i_local; lam_p2_local];
        else
            q_pair = group.q_pair;
            rout_fine_other2 = [q_pair(2)+pair_cache.rout_base_f; group.refine_canon{2}];
            R2 = -lapSLPmat(q_pair(1)+rbase_in_c,rout_fine_other2)*lam_i_local;

            rout_fine_other1 = [q_pair(1)+pair_cache.rout_base_f; group.refine_canon{1}];
            R1 = -lapSLPmat(q_pair(2)+rbase_in_c,rout_fine_other1)*lam_p2_local;

            rhs_pair = [R1; R2];
        end

        pair_mapped = group.Upf*rhs_pair;
        beta_tot_nonp_local = group.Ypf*pair_mapped;
        if isfield(opt,'cmap') && opt.cmap
            tau_peanut_nonp_local = group.Cmap*rhs_pair;
        else
            tau_peanut_nonp_local = group.YC*(group.DC*beta_tot_nonp_local);
        end

        rimage_i = pair.rimage_i;
        rimage_p2 = pair.rimage_j;
        im_i = numel(group.rimage_canon{1});
        im_p2 = numel(group.rimage_canon{2});
    else
        if precomp
            rhs_pair = [lam_i; lam_p2];
        else
            rout_fine_other2 = getFineOther(opt.a_f,opt.N_f,refine,q,i,p2,opt.rad);
            R2 = -lapSLPmat(rbase_in_c+q(i),rout_fine_other2)*lam_i;

            rout_fine_other1 = getFineOther(opt.a_f,opt.N_f,refine,q,p2,i,opt.rad);
            R1 = -lapSLPmat(rbase_in_c+q(p2),rout_fine_other1)*lam_p2;

            rhs_pair = [R1; R2];
        end

        pair_mapped = Upf{i,p2}*rhs_pair;
        beta_tot_nonp_local = Ypf{i,p2}*pair_mapped;

        if isfield(opt,'cmap') && opt.cmap
            tau_peanut_nonp_local = Cmap{i,p2}*rhs_pair;
        else
            tau_peanut_nonp_local = YC_all{i,p2}*(DC_all{i,p2}*beta_tot_nonp_local);
        end

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
    beta_i_local = projectChargeMode(beta_i_nonp_local,project_charge);
    beta_p2_local = projectChargeMode(beta_p2_nonp_local,project_charge);
    tau_peanut_local = [projectChargeMode(tau_peanut_nonp_local(1:N_c),project_charge); ...
        projectChargeMode(tau_peanut_nonp_local(N_c+1:2*N_c),project_charge)];

    if use_pair_cache
        beta_i_fine = rotateUniformCircleData(beta_i_local(1:N_f),conj(rot));
        beta_p2_fine = rotateUniformCircleData(beta_p2_local(1:N_f),conj(rot));
        beta_i_fine_nonp = rotateUniformCircleData(beta_i_nonp_local(1:N_f),conj(rot));
        beta_p2_fine_nonp = rotateUniformCircleData(beta_p2_nonp_local(1:N_f),conj(rot));

        tau_peanut_i = rotateUniformCircleData(tau_peanut_local(1:N_c),conj(rot));
        tau_peanut_p2 = rotateUniformCircleData(tau_peanut_local(N_c+1:2*N_c),conj(rot));
        tau_peanut_i_nonp = rotateUniformCircleData(tau_peanut_nonp_local(1:N_c),conj(rot));
        tau_peanut_p2_nonp = rotateUniformCircleData(tau_peanut_nonp_local(N_c+1:2*N_c),conj(rot));
    else
        beta_i_fine = beta_i_local(1:N_f);
        beta_p2_fine = beta_p2_local(1:N_f);
        beta_i_fine_nonp = beta_i_nonp_local(1:N_f);
        beta_p2_fine_nonp = beta_p2_nonp_local(1:N_f);

        tau_peanut_i = tau_peanut_local(1:N_c);
        tau_peanut_p2 = tau_peanut_local(N_c+1:2*N_c);
        tau_peanut_i_nonp = tau_peanut_nonp_local(1:N_c);
        tau_peanut_p2_nonp = tau_peanut_nonp_local(N_c+1:2*N_c);
    end

    tau_peanut_actual = [tau_peanut_i; tau_peanut_p2];
    tau_peanut_nonp_actual = [tau_peanut_i_nonp; tau_peanut_p2_nonp];
    beta_tot_actual = [beta_i_fine; beta_i_local(N_f+1:end); beta_p2_fine; beta_p2_local(N_f+1:end)];

    lam_c(idx_i) = lam_c(idx_i) + tau_peanut_i;
    lam_c(idx_p2) = lam_c(idx_p2) + tau_peanut_p2;
    lam_c_nonp(idx_i) = lam_c_nonp(idx_i) + tau_peanut_i_nonp;
    lam_c_nonp(idx_p2) = lam_c_nonp(idx_p2) + tau_peanut_p2_nonp;

    lam_f{i} = lam_f{i} + beta_i_fine;
    lam_f{p2} = lam_f{p2} + beta_p2_fine;
    lam_e_chunks{i}{end+1,1} = beta_i_local(N_f+1:end);
    lam_e_chunks{p2}{end+1,1} = beta_p2_local(N_f+1:end);
    lam_f_nonp{i} = lam_f_nonp{i} + beta_i_fine_nonp;
    lam_f_nonp{p2} = lam_f_nonp{p2} + beta_p2_fine_nonp;
    lam_e_nonp_chunks{i}{end+1,1} = beta_i_nonp_local(N_f+1:end);
    lam_e_nonp_chunks{p2}{end+1,1} = beta_p2_nonp_local(N_f+1:end);

    block_i = (i-1)*N_check+1:i*N_check;
    block_p2 = (p2-1)*N_check+1:p2*N_check;
    rout_pair = [rcheck_out(block_i); rcheck_out(block_p2)];

    rin_pair_f = [q(i)+rbase_in_f; rimage_i; q(p2)+rbase_in_f; rimage_p2];
    rin_pair_c = [q(i)+rbase_in_c; q(p2)+rbase_in_c];

    u_fine = lapSLPfield(rin_pair_f,rout_pair,beta_tot_actual,false);
    u_peanut = lapSLPfield(rin_pair_c,rout_pair,tau_peanut_actual,false);

    pair_idx = [block_i block_p2]';
    u_corr(pair_idx) = u_corr(pair_idx) + u_fine - u_peanut;
end

for k = 1:P
    if ~isempty(lam_e_chunks{k})
        lam_e{k} = vertcat(lam_e_chunks{k}{:});
    end
    if ~isempty(lam_e_nonp_chunks{k})
        lam_e_nonp{k} = vertcat(lam_e_nonp_chunks{k}{:});
    end
end

end

function lam_out = projectChargeMode(lam_in,project_charge)
%PROJECTCHARGEMODE Apply scalar charge projection with Kq = ones. Only
%needed for the elastance problem.

if ~project_charge || isempty(lam_in)
    lam_out = lam_in;
    return
end

n = numel(lam_in);
lam_out = lam_in - (sum(lam_in)/n);

end

function pair_cache = get_pair_cache(geom,basis)
if isfield(basis,'pair_cache') && ~isempty(basis.pair_cache)
    pair_cache = basis.pair_cache;
elseif isfield(geom,'pair_cache') && ~isempty(geom.pair_cache)
    pair_cache = geom.pair_cache;
else
    pair_cache = [];
end
end
