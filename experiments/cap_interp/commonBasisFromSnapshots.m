function basis = commonBasisFromSnapshots(dA_train,r_max)
%COMMONBASISFROMSNAPSHOTS Common left/right spatial bases for gap-dependent part.
%
% basis = commonBasisFromSnapshots(dA_train,r_max) forms
%   X_L = [dA(:,:,1), ..., dA(:,:,q)]           (n x n*q)
%   X_R = [dA(:,:,1)', ..., dA(:,:,q)']         (n x n*q)
% from the training snapshots dA_train = A(s_j)-A_ref (n x n x q), and
% returns their economical SVDs truncated to r_max columns.
%
% basis.U, basis.sigma_L   - left common basis (n x r_max) and its
%                            singular values (r_max x 1)
% basis.V, basis.sigma_R   - right common basis (n x r_max) and its
%                            singular values (r_max x 1)
%
% U(:,1:r_L) and V(:,1:r_R) for any r_L,r_R <= r_max are themselves valid
% truncated common bases (nested, since they come from one SVD).
%
% Anna Broms, Sep 2026

production_basis = buildLaplaceCmapSpatialBasis(dA_train,r_max);
basis.U = production_basis.U;
basis.V = production_basis.V;
basis.sigma_L = production_basis.sigma_left;
basis.sigma_R = production_basis.sigma_right;

end
