% t = clustered_points(N, Ns[, sigma, A])
%
% Compute N+2Ns points on the interval [-1,1]:
% - N points are equispaced on the interval (excluding endpoints)
% - Ns points cluster towards the left and Ns points towards the right,
%   on the intervals [-1..-1+A] and [1-A..1] respectively.
%
% Input:
% - N: number of equispaced points
% - Ns: number of clustering points
% - sigma: distribution parameter (default: 4)
% - A: scaling of clusters (default: 1/2)
%
% Output:
% - t: the points
function t = clustered_points(N, Ns, sigma, A)
    if nargin < 3
        sigma = 4;
    end
    if nargin < 4
        A = 1/2;
    end
    t1 = linspace(-1, 1, 2*N+1)';
    t1 = t1(2:2:end-1);

    u = sqrt(1:Ns)' - sqrt(Ns);
    v = A*exp(sigma*u);
    t2 = [-1+v; 1-v];

    t = sort([t1; t2]);
end
