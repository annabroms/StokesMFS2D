function basis = buildLaplaceCmapSpatialBasis(dC_snapshots,r_max)
%BUILDLAPLACECMAPSPATIALBASIS Common left/right bases for Cmap increments.

[n,~,q] = size(dC_snapshots);
X_left = reshape(dC_snapshots,n,n*q);
X_right = reshape(permute(dC_snapshots,[2 1 3]),n,n*q);
r_max = min([r_max,n,n*q]);

[U,S_left,~] = svd(X_left,'econ');
[V,S_right,~] = svd(X_right,'econ');
basis = struct();
basis.U = U(:,1:r_max);
basis.V = V(:,1:r_max);
sigma_left = diag(S_left);
sigma_right = diag(S_right);
basis.sigma_left = sigma_left(1:r_max);
basis.sigma_right = sigma_right(1:r_max);
end
