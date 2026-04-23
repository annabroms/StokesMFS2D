function res = matvec_lap_peanut_big_sparse(tau,geom,basis)
%MATVEC_LAP_PEANUT_BIG_SPARSE Laplace peanut matvec using sparse pair maps.

validateBigSparseMatvecOptions(geom,basis);

rvec_in = geom.rvec_in;
opt = geom.opt;
pairs = geom.pairs;
rcheck = geom.rcheck;

P = opt.P;
N_c = opt.N_c;
N_large = round(N_c*opt.a_c);

[lam_c,lam_self,u_corr,pair_qv_nonp] = ...
    transform_lap_peanut_big_sparse(tau,geom,basis);

res = lapSLPfield(rvec_in,rcheck,lam_c,opt.use_fmm);
res = res + u_corr;

if isfield(opt,'project_charge') && logical(opt.project_charge)
    has_neigh = false(1,P);
    has_neigh(pairs(:)) = true;
    if any(has_neigh)
        % Add the same pair-induced charge correction across each touched
        % body block in one batch instead of a MATLAB loop.
        res_blocks = reshape(res,N_large,P);
        res_blocks(:,has_neigh) = bsxfun(@plus, ...
            res_blocks(:,has_neigh), ...
            reshape(pair_qv_nonp(has_neigh),1,[]));
        res = res_blocks(:);

        % % this is the same as
        % has_neigh = sort(unique(pairs(:)));
        % for ii = 1:numel(has_neigh)
        %     k = has_neigh(ii);
        %     block = (k-1)*N_large+1:k*N_large;
        %     res(block) = res(block) + pair_qv_nonp(k);
        % end
    end
end

% Correct all diagonal identity blocks in one batch instead of a more
% expensive loop over the P particles.
lam_self_blocks = reshape(lam_self,N_c,P);
uii_blocks = basis.Nii*lam_self_blocks;
res = res - reshape(uii_blocks(1:N_large,:),[],1);

res = res + tau;

end

function validateBigSparseMatvecOptions(geom,basis)
opt = geom.opt;
if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('matvec_lap_peanut_big_sparse:MissingBigSparse', ...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('matvec_lap_peanut_big_sparse:UnsupportedGrid', ...
        'The Laplace big-sparse matvec is only supported on the solve grid.');
end
if ~logical(getOptField(opt,'cmap',false))
    error('matvec_lap_peanut_big_sparse:UnsupportedOption', ...
        'opt.use_big_sparse=1 requires opt.cmap=1.');
end

end
