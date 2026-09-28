function r = numericalRank(sigma,tol)
%NUMERICALRANK Number of singular values above tol*sigma(1).
%
% r = numericalRank(sigma,tol), sigma a descending singular-value vector.
% tol may be a vector, in which case r is returned elementwise.
%
% Anna Broms, Sep 2026

sigma = sigma(:);
if isempty(sigma) || sigma(1)==0
    r = zeros(size(tol));
    return
end
r = arrayfun(@(t) sum(sigma > t*sigma(1)),tol);

end
