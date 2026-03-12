function res = matvec_lap_peanut_enhanced(tau,geom,basis)
%MATVEC_LAP_PEANUT_ENHANCED Matrix-vector action for Laplace peanut preconditioner.
%
% Syntax:
%   res = matvec_lap_peanut_enhanced(tau,geom,basis)
%
% See also: solve_cap_peanut, transform_lap_peanut.
%
% Anna Broms, Mar 2026

rbase_in_c = geom.rbase_in_c;
rvec_in = geom.rvec_in;
rvec_out = geom.rvec_out;
opt = geom.opt;
q = geom.q;
pairs = geom.pairs;

if isfield(geom,'rcheck') && ~isempty(geom.rcheck)
    rcheck = geom.rcheck;
else
    rcheck = rvec_out;
end

P = numel(q);
N_large = length(rvec_out)/P;

%%  Get source strengths from data on boundary
[lam_c,lam_self,~,~,u_corr,~,~,lam_f_nonp,lam_e_nonp] = ...
    transform_lap_peanut(tau,geom,basis);

%% Evaluate field
res = lapSLPField(rvec_in,rcheck,lam_c,opt.use_fmm);
res = res + u_corr;


if isequal(rcheck,rvec_out)

    %% For elastance: add Lr blocks on diagonal from fine sources
    if isfield(opt,'project_charge') && logical(opt.project_charge)
        has_neigh = sort(unique(pairs(:)));
        for ii = 1:numel(has_neigh)
            k = has_neigh(ii);
            block = (k-1)*N_large+1:k*N_large;
            lr_pair = sum(lam_f_nonp{k}) + sum(lam_e_nonp{k});
            res(block) = res(block) + lr_pair;
        end
    end
    
    %% Correct identity diagonal blocks
    rout = rvec_out(1:N_large)-q(1);
    Nii = lapSLPmat(rbase_in_c,rout);

    for i = 1:P
        idx = (i-1)*opt.N_c+1:i*opt.N_c;
        block = (i-1)*N_large+1:i*N_large;
        uii = Nii*lam_self(idx);
        res(block) = res(block)-uii;
    end

    res = res + tau;
end

end
