function [lam_c,lam_self,lam_f,lam_e,u_corr,pair_qv_nonp,lam_c_nonp,lam_self_nonp,lam_f_nonp,lam_e_nonp] = transform_lap_peanut(tau,geom,basis)
%TRANSFORM_LAP_PEANUT Map coarse boundary data to interior compressed/fine
%Laplace sources. Helper function to matvec_lap_peanut_enhanced
%
% Syntax:
%   [lam_c,lam_self,lam_f,lam_e,u_corr] = transform_lap_peanut(tau,geom,basis)
%   [lam_c,lam_self,lam_f,lam_e,u_corr,pair_qv_nonp,lam_c_nonp, ...
%       lam_self_nonp,lam_f_nonp,lam_e_nonp] = transform_lap_peanut(tau,geom,basis)
%
% Outputs:
%   lam_*      - Projected source strengths.
%   pair_qv_nonp
%              - Per-body sums of the non-projected pair-induced fine/image
%                sources. For cmap=1 this is recovered via basis.Cmap_QV.
%   lam_*_nonp - Non-projected source strengths (used for Lr/voltage terms).
%                Fine/image outputs are returned empty when explicit pair
%                sources are not needed.
%
% See also: getPairBasisLaplace, getPeanutBlockLaplace, ...
%   matvec_lap_peanut_enhanced.
%
% Anna Broms, Mar 2026

rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
rimage_vec = geom.rimage_vec;
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
Cmap_QV = basis.Cmap_QV;
pair_cache = basis.pair_cache;
use_pair_cache = pair_cache.enabled;
use_cmap = isfield(opt,'cmap') && opt.cmap;
get_bndry_field = logical(getOptField(opt,'get_bndry_field',true));
need_explicit_pair_sources = ~use_cmap || get_bndry_field;

P = numel(q);
N_c = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P;
N_check = length(rcheck_out)/P;
use_ucorr_colloc = use_pair_cache && isequal(rcheck_out,rvec_out);
if isfield(opt,'project_charge') && ~isempty(opt.project_charge)
    project_charge = logical(opt.project_charge);
else
    project_charge = false;
end

pair_qv_nonp = zeros(P,1);
if need_explicit_pair_sources
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
else
    lam_f = [];
    lam_e = [];
    lam_f_nonp = [];
    lam_e_nonp = [];
end

u_corr = zeros(P*N_check,1);

% Phase 1: one-body coarse mapping.
% for i = 1:P
%     block = (i-1)*N_large+1:i*N_large;
%     tau_i = tau(block);

%     lam_i_nonp = Y{1}*(U{1}*tau_i);
%     lam_i = projectChargeMode(lam_i_nonp,project_charge);

%     idx = (i-1)*N_c+1:i*N_c;
%     lam_c(idx) = lam_i;
%     lam_c_nonp(idx) = lam_i_nonp;
% end

% Faster alternative: All one-body blocks are identical, so batch them as
% particle columns while keeping the U{1} then Y{1} ordering for a backward stable apply.
tau_blocks = reshape(tau(1:P*N_large),N_large,P);
lam_blocks_nonp = Y{1}*(U{1}*tau_blocks);
if project_charge
    lam_blocks = lam_blocks_nonp - mean(lam_blocks_nonp,1);
else
    lam_blocks = lam_blocks_nonp;
end
lam_c = reshape(lam_blocks,[],1);
lam_c_nonp = reshape(lam_blocks_nonp,[],1);

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

    if use_pair_cache % if distances are shared for more than one pair
        meta = pair_cache.meta(row);
        group = pair_cache.groups(meta.group_id);
        
        % interpolate data to canonical reference frame
        rhs_pair = rotateUniformCircleData([lam_i lam_p2],[],meta.phase_c);
        rhs_pair = rhs_pair(:);

        % recover fine sources
        if need_explicit_pair_sources
            pair_mapped = group.Upf*rhs_pair;
            beta_tot_nonp_local = group.Ypf*pair_mapped;
        else
            beta_tot_nonp_local = [];
        end

        % use coarse-to-coarse map without going via fine sources
        if use_cmap
            tau_peanut_nonp_local = group.Cmap*rhs_pair;
            pair_qv_local = group.Cmap_QV*rhs_pair;
        else
            tau_peanut_nonp_local = group.YC*(group.DC*beta_tot_nonp_local);
            pair_qv_local = [];
        end
    else
        rhs_pair = [lam_i; lam_p2];
        pair_mapped = Upf{i,p2}*rhs_pair;
        if need_explicit_pair_sources
            beta_tot_nonp_local = Ypf{i,p2}*pair_mapped;
        else
            beta_tot_nonp_local = [];
        end

        if use_cmap
            tau_peanut_nonp_local = Cmap{i,p2}*rhs_pair;
            pair_qv_local = Cmap_QV{i,p2}*rhs_pair;
        else
            tau_peanut_nonp_local = YC_all{i,p2}*(DC_all{i,p2}*beta_tot_nonp_local);
            pair_qv_local = [];
        end
    end

    tau_peanut_nonp_pair = reshape(tau_peanut_nonp_local,N_c,2);
   % tau_peanut_pair = [projectChargeMode(tau_peanut_nonp_pair(:,1),project_charge) ...
      %  projectChargeMode(tau_peanut_nonp_pair(:,2),project_charge)];

    if use_pair_cache 
        % tau_pair_rot = rotateUniformCircleData([tau_peanut_pair tau_peanut_nonp_pair],[],meta.phase_c_inv);
        % tau_peanut_pair = tau_pair_rot(:,1:2);
        % tau_peanut_nonp_pair = tau_pair_rot(:,3:4);
        tau_pair_rot = rotateUniformCircleData(tau_peanut_nonp_pair,[],meta.phase_c_inv);
        tau_peanut_nonp_pair = tau_pair_rot;
    end
    tau_peanut_pair = [projectChargeMode(tau_peanut_nonp_pair(:,1),project_charge) ...
        projectChargeMode(tau_peanut_nonp_pair(:,2),project_charge)];

    lam_c(idx_i) = lam_c(idx_i) + tau_peanut_pair(:,1);
    lam_c(idx_p2) = lam_c(idx_p2) + tau_peanut_pair(:,2);
    lam_c_nonp(idx_i) = lam_c_nonp(idx_i) + tau_peanut_nonp_pair(:,1);
    lam_c_nonp(idx_p2) = lam_c_nonp(idx_p2) + tau_peanut_nonp_pair(:,2);
    if use_cmap
        pair_qv_nonp(i) = pair_qv_nonp(i) + pair_qv_local(1);
        pair_qv_nonp(p2) = pair_qv_nonp(p2) + pair_qv_local(2);
    end


    %% Determine correction field for the pair
    block_i = (i-1)*N_check+1:i*N_check;
    block_p2 = (p2-1)*N_check+1:p2*N_check;
    rout_pair = [rcheck_out(block_i); rcheck_out(block_p2)];
    rin_pair_c = [q(i)+rbase_in_c; q(p2)+rbase_in_c];
    use_dense_u_pair = use_ucorr_colloc && use_pair_cache && ~need_explicit_pair_sources ...
        && isfield(meta,'Ucross_colloc_actual') && ~isempty(meta.Ucross_colloc_actual) ...
        && isfield(meta,'Ec_colloc_actual') && ~isempty(meta.Ec_colloc_actual);    

    if use_dense_u_pair
        % use stored dense matrices
        rhs_actual = [lam_i; lam_p2];
        u_fine = meta.Ucross_colloc_actual*rhs_actual;

        if project_charge
            u_fine = u_fine - meta.Lr_colloc_actual*pair_qv_local;
        end

        u_peanut = meta.Ec_colloc_actual*tau_peanut_pair(:);
        u_pair = u_fine - u_peanut;
    elseif need_explicit_pair_sources 
        % use fine sources
        if use_pair_cache
            rimage_i = meta.mid + meta.rot*group.rimage_canon{1};
            rimage_p2 = meta.mid + meta.rot*group.rimage_canon{2};
        else
            rimage_i = rimage_vec{i,p2};
            rimage_p2 = rimage_vec{p2,i};
        end
        im_i = numel(rimage_i);
        im_p2 = numel(rimage_p2);
        s_i = 1:N_f;
        e_i = N_f+1:N_f+im_i;
        s_p2 = N_f+im_i+1:2*N_f+im_i;
        e_p2 = 2*N_f+im_i+1:2*N_f+im_i+im_p2;
        beta_i_nonp_local = [beta_tot_nonp_local(s_i); beta_tot_nonp_local(e_i)];
        beta_p2_nonp_local = [beta_tot_nonp_local(s_p2); beta_tot_nonp_local(e_p2)];
        
        beta_i_local = projectChargeMode(beta_i_nonp_local,project_charge);
        beta_p2_local = projectChargeMode(beta_p2_nonp_local,project_charge);

        if ~use_cmap
            pair_qv_local = [sum(beta_i_nonp_local); sum(beta_p2_nonp_local)];
            pair_qv_nonp(i) = pair_qv_nonp(i) + pair_qv_local(1);
            pair_qv_nonp(p2) = pair_qv_nonp(p2) + pair_qv_local(2);
        end

        beta_fine_pair = [beta_i_local(1:N_f) beta_p2_local(1:N_f) ...
            beta_i_nonp_local(1:N_f) beta_p2_nonp_local(1:N_f)];
        if use_pair_cache
            beta_fine_pair = rotateUniformCircleData(beta_fine_pair,[],meta.phase_f_inv);
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

        beta_tot_actual = [beta_i_fine; beta_i_local(N_f+1:end); ...
            beta_p2_fine; beta_p2_local(N_f+1:end)];
        rin_pair_f = [q(i)+rbase_in_f; rimage_i; q(p2)+rbase_in_f; rimage_p2];
        u_fine = lapSLPfield(rin_pair_f,rout_pair,beta_tot_actual,false);
        u_peanut = lapSLPfield(rin_pair_c,rout_pair,tau_peanut_pair(:),false);
        u_pair = u_fine - u_peanut;

    else
        % use direct summation to recreate rhs
        u_fine_i = -lapSLPfield(q(p2)+rbase_in_c,rcheck_out(block_i),lam_p2,false);
        u_fine_p2 = -lapSLPfield(q(i)+rbase_in_c,rcheck_out(block_p2),lam_i,false);

        % To match rhs in pair-problem, ones-matrix must be subtracted as
        % it is later added again. 
        if project_charge
            u_fine_i = u_fine_i - pair_qv_local(1);
            u_fine_p2 = u_fine_p2 - pair_qv_local(2);
        end
        u_fine = [u_fine_i; u_fine_p2];
        u_peanut = lapSLPfield(rin_pair_c,rout_pair,tau_peanut_pair(:),false);
        u_pair = u_fine - u_peanut;
    end

    pair_idx = [block_i block_p2]';
    u_corr(pair_idx) = u_corr(pair_idx) + u_pair;
end


%% Store sources
if need_explicit_pair_sources
    for k = 1:P
        if ~isempty(lam_e_chunks{k})
            lam_e{k} = vertcat(lam_e_chunks{k}{:});
        end
        if ~isempty(lam_e_nonp_chunks{k})
            lam_e_nonp{k} = vertcat(lam_e_nonp_chunks{k}{:});
        end
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

