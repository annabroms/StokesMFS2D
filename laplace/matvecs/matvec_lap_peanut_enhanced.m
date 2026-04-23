function res = matvec_lap_peanut_enhanced(tau,geom,basis)
%MATVEC_LAP_PEANUT_ENHANCED Matrix-vector action for Laplace peanut preconditioner.
%
% Syntax:
%   res = matvec_lap_peanut_enhanced(tau,geom,basis)
%
% See also: solve_cap_peanut, transform_lap_peanut.
%
% Anna Broms, Mar 2026

rvec_in = geom.rvec_in;
rvec_out = geom.rvec_out;
opt = geom.opt;
pairs = geom.pairs;

if isfield(geom,'rcheck') && ~isempty(geom.rcheck)
    rcheck = geom.rcheck;
else
    rcheck = rvec_out;
end

P = opt.P;
N_c = opt.N_c;
N_large = round(N_c*opt.a_c);

%%  Get source strengths from data on boundary
[lam_c,lam_self,~,~,u_corr,pair_qv_nonp] = ...
    transform_lap_peanut(tau,geom,basis);

%% Evaluate field
res = lapSLPfield(rvec_in,rcheck,lam_c,opt.use_fmm);
res = res + u_corr;


if isequal(rcheck,rvec_out)

    %% For elastance: add Lr blocks on diagonal from fine sources
    if isfield(opt,'project_charge') && logical(opt.project_charge)
        has_neigh = false(1,P);
        has_neigh(pairs(:)) = true;
        if any(has_neigh)
            % Add the same pair-induced charge correction across each
            % touched body block in one batch instead of a MATLAB loop.
            res_blocks = reshape(res,N_large,P);
            res_blocks(:,has_neigh) = bsxfun(@plus, ...
                res_blocks(:,has_neigh), ...
                reshape(pair_qv_nonp(has_neigh),1,[]));
            res = res_blocks(:);
        end
    end
    
    % Correct all diagonal identity blocks in one batch instead of a more
    % expensive loop over the P particles.

    lam_self_blocks = reshape(lam_self,N_c,P);
    uii_blocks = basis.Nii*lam_self_blocks;
    res = res - reshape(uii_blocks(1:N_large,:),[],1);

    res = res + tau;
end

end
