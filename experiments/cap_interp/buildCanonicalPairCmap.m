function [C,info] = buildCanonicalPairCmap(delta,R,opt,grids)
%BUILDCANONICALPAIRCMAP Exact canonical aligned pair-correction matrix C(delta).
%
% [C,info] = buildCanonicalPairCmap(delta,R,opt,grids) builds the same
% fixed-size coarse-to-coarse pair map Cmap{1,2} that is inserted into the
% global coarse Laplace system for a canonical, aligned pair of two equal
% discs of radius R with gap delta, centers on the real axis at 0 and
% 2R+delta. This is a thin wrapper around the production routines
% getEnhancedGrid and getPairBasisLaplace (opt.cmap=1), reusing exactly the
% construction used by sep24_capacitance_cmap_delta_smoothness.m and
% sep25_capacitance_cmap_rank_truncation.m. No approximation or
% interpolation is done here; every call is an exact build.
%
% grids is a struct with fields rbase_in_c, rbase_in_f, rout_base_f,
% rbase_out_c (fixed discretisation grids, independent of delta), as built
% in getLaplace2Dparams-derived scripts.
%
% info fields include effective_dof_fine/effective_dof_peanut, the traces
% of the two Tikhonov influence matrices. The legacy rank_fine/rank_peanut
% fields report what a 1e-14 TSVD diagnostic would retain; they do not
% control a Tikhonov build. n_image is the active ellipse-source count,
% t_build times only the production Cmap build, and n is the matrix size.
%
% Anna Broms, Sep 2026

t0 = tic;
q = [0;2*R+delta];
[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
assert(isequal(pairs,[1 2]),'buildCanonicalPairCmap:BadPairs', ...
    'Expected a single close pair (1,2) for this canonical aligned setup.');

snapshot = buildCanonicalLaplacePairMap(delta,opt,grids);
C = snapshot.Cmap;
t_build = toc(t0);

info = struct();
info.n = size(C,1);

% Recreate the two production matrices only for diagnostics. These SVDs do
% not affect C and are deliberately excluded from t_build above.
rin_pair_f = [q(1)+grids.rbase_in_f; rimage_vec{1,2}; ...
    q(2)+grids.rbase_in_f; rimage_vec{2,1}];
rout_f = [q(1)+grids.rout_base_f; refine{1,2}; ...
    q(2)+grids.rout_base_f; refine{2,1}];
s_fine = svd(lapSLPmat(rin_pair_f,rout_f));

% getPeanutBlockLaplace's rank is not returned when opt.cmap=1, so replicate
% the same matrix it factorises purely to read off its numerical rank (no
% effect on the production Cmap computation), exactly as in sep24.
rin_pair_c = [q(1)+grids.rbase_in_c; q(2)+grids.rbase_in_c];
rout_peanut = createPeanut(q(1),q(2),opt.N_peanut,0,R);
s_peanut = svd(lapSLPmat(rin_pair_c,rout_peanut));
svd_tol = 1e-14; % matches the hardcoded tolerance in both blocks (sep24)
info.rank_fine = sum(s_fine > max(s_fine)*svd_tol);
info.rank_peanut = sum(s_peanut > max(s_peanut)*svd_tol);
info.effective_dof_fine = regularized_dof(s_fine,opt,svd_tol);
info.effective_dof_peanut = regularized_dof(s_peanut,opt,svd_tol);
info.n_image = [numel(rimage_vec{1,2}) numel(rimage_vec{2,1})];
info.t_build = t_build;

end

function dof = regularized_dof(s,opt,svd_tol)
if isempty(s) || max(s)==0
    dof = 0;
elseif logical(getOptField(opt,'use_tikhonov',false))
    eta = getOptField(opt,'tikhonov_tol',svd_tol);
    sigma_scaled = s/max(s);
    dof = sum(sigma_scaled.^2./(sigma_scaled.^2+eta^2));
else
    dof = sum(s > max(s)*svd_tol);
end
end
