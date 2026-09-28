function Fq = evalMatrixChebBary(s_query,s_nodes,w,F_nodes)
%EVALMATRIXCHEBBARY Barycentric Chebyshev interpolation of matrix snapshots.
%
% Fq = evalMatrixChebBary(s_query,s_nodes,w,F_nodes) evaluates the unique
% degree-(q-1) polynomial interpolant (applied entrywise, i.e. one scalar
% interpolant per matrix entry, all sharing the same nodes/weights) through
% the snapshots F_nodes(:,:,j) at s_nodes(j), at each point in s_query.
%
% s_nodes, w  - q x 1, from chebLobattoNodes/chebBarycentricWeights (same q).
% F_nodes     - n1 x n2 x q array of matrix snapshots.
% s_query     - m x 1 (or scalar) evaluation points.
%
% Fq is n1 x n2 x m. Uses the 2nd-kind barycentric formula (Berrut &
% Trefethen 2004), which is numerically stable and reproduces a training
% snapshot exactly (to round-off) when s_query coincides with a node.
%
% Anna Broms, Sep 2026

s_query = s_query(:);
q = numel(s_nodes);
[n1,n2,q_check] = size(F_nodes);
assert(q_check==q,'evalMatrixChebBary:SizeMismatch', ...
    'F_nodes must have one snapshot per node in s_nodes.');

m = numel(s_query);
F_flat = reshape(F_nodes,n1*n2,q);
Fq_flat = zeros(n1*n2,m);

exact_tol = 1e-12*max(1,max(abs(s_nodes)));
for k = 1:m
    d = s_query(k) - s_nodes;
    [dmin,j_exact] = min(abs(d));
    if dmin < exact_tol
        Fq_flat(:,k) = F_flat(:,j_exact);
        continue
    end
    c = w./d;
    Fq_flat(:,k) = (F_flat*c)/sum(c);
end

Fq = reshape(Fq_flat,n1,n2,m);

end
