function res = matvec_lap_2B_enhanced(tau,geom,basis,rcheck)
%MATVEC_LAP_2B_ENHANCED Matrix-vector action for Laplace 2B preconditioner.
%
% Syntax:
%   res = matvec_lap_2B_enhanced(tau,geom,basis,rcheck)
%
% See also: solve_cap_2B, solve_elast_2B, getPairTransformationLaplace.
%
% Anna Broms, Mar 2026

rbase_in_c = geom.rbase_in_c;
rvec_out = geom.rvec_out;
q = geom.q;
opt = geom.opt;
if nargin < 4 || isempty(rcheck)
    rcheck = rvec_out;
end

P = numel(q);
N_large = length(rvec_out)/P;

[rvec_in,coarse_ind,lam_all,lam_c,~,~,~,~,lam_f_nonp,lam_e_nonp] = ...
    getPairTransformationLaplace(tau,geom,basis);

res = lapSLPField(rvec_in,rcheck,lam_all,opt.use_fmm);

if isequal(rcheck,rvec_out)
    if isfield(opt,'project_charge') && logical(opt.project_charge)
        % Add the Lr action for pair-induced source components (mobility-style).
        has_neigh = sort(unique(geom.pairs(:)));
        for ii = 1:numel(has_neigh)
            k = has_neigh(ii);
            block = (k-1)*N_large+1:k*N_large;
            lr_pair = sum(lam_f_nonp{k}) + sum(lam_e_nonp{k});
            res(block) = res(block) + lr_pair;
        end
    end

    rout = rvec_out(1:N_large)-q(1);
    Nii = lapSLPmat(rbase_in_c,rout);

    for i = 1:P
        idx = coarse_ind{i};
        uii = Nii*lam_c(idx);
        block = (i-1)*N_large+1:i*N_large;
        res(block) = res(block)-uii;
    end

    res = res+tau;
end

end
