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
    beta_tot_nonp = Ypf{i,p2}*pair_mapped;

    rimage_i = rimage_vec{i,p2};
    rimage_p2 = rimage_vec{p2,i};
    im_i = length(rimage_i);
    im_p2 = length(rimage_p2);

    s_i = 1:N_f;
    e_i = N_f+1:N_f+im_i;
    s_p2 = N_f+im_i+1:2*N_f+im_i;
    e_p2 = 2*N_f+im_i+1:2*N_f+im_i+im_p2;

    beta_i_nonp = [beta_tot_nonp(s_i); beta_tot_nonp(e_i)];
    beta_p2_nonp = [beta_tot_nonp(s_p2); beta_tot_nonp(e_p2)];
    beta_i = projectChargeMode(beta_i_nonp,project_charge);
    beta_p2 = projectChargeMode(beta_p2_nonp,project_charge);
    beta_tot = [beta_i(1:N_f); beta_i(N_f+1:end); beta_p2(1:N_f); beta_p2(N_f+1:end)];

    if isfield(opt,'cmap') && opt.cmap
        tau_peanut_nonp = Cmap{i,p2}*rhs_pair;
    else
        tau_peanut_nonp = YC_all{i,p2}*(DC_all{i,p2}*beta_tot_nonp);
    end
    tau_peanut = [projectChargeMode(tau_peanut_nonp(1:N_c),project_charge); ...
        projectChargeMode(tau_peanut_nonp(N_c+1:2*N_c),project_charge)];

    lam_c(idx_i) = lam_c(idx_i) + tau_peanut(1:N_c);
    lam_c(idx_p2) = lam_c(idx_p2) + tau_peanut(N_c+1:2*N_c);
    lam_c_nonp(idx_i) = lam_c_nonp(idx_i) + tau_peanut_nonp(1:N_c);
    lam_c_nonp(idx_p2) = lam_c_nonp(idx_p2) + tau_peanut_nonp(N_c+1:2*N_c);

    lam_f{i} = lam_f{i} + beta_i(1:N_f);
    lam_f{p2} = lam_f{p2} + beta_p2(1:N_f);
    lam_e_chunks{i}{end+1,1} = beta_i(N_f+1:end);
    lam_e_chunks{p2}{end+1,1} = beta_p2(N_f+1:end);
    lam_f_nonp{i} = lam_f_nonp{i} + beta_i_nonp(1:N_f);
    lam_f_nonp{p2} = lam_f_nonp{p2} + beta_p2_nonp(1:N_f);
    lam_e_nonp_chunks{i}{end+1,1} = beta_i_nonp(N_f+1:end);
    lam_e_nonp_chunks{p2}{end+1,1} = beta_p2_nonp(N_f+1:end);

    block_i = (i-1)*N_check+1:i*N_check;
    block_p2 = (p2-1)*N_check+1:p2*N_check;
    rout_pair = [rcheck_out(block_i); rcheck_out(block_p2)];

    rin_pair_f = [q(i)+rbase_in_f; rimage_i; q(p2)+rbase_in_f; rimage_p2];
    rin_pair_c = [q(i)+rbase_in_c; q(p2)+rbase_in_c];

    u_fine = lapSLPField(rin_pair_f,rout_pair,beta_tot,false);
    u_peanut = lapSLPField(rin_pair_c,rout_pair,tau_peanut,false);

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
