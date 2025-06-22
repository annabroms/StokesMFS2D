function t = oversample_points(x, L, a, b)
%function t = oversample_points(x, L, a, b)
%
%   Return a vector of points based on the given points x, but with an
%   additional L-1 points between each two elements of x. The elements of
%   x are included in the output. It is assumed that x is sorted.
%
%   Optionally, an interval [a,b] is specified, and L-1 points are included
%   in the interval [a,x(1)] and [x(end),b].
%
%   The resulting vector is approximately L times larger than x.

t = [x(1)];
for i = 1:length(x)-1
    ti = linspace(x(i), x(i+1), L+1);
    ti = ti(2:end);
    t = [t ti];
end

if nargin == 4
    ti = linspace(a, x(1), L+1);
    ti = ti(2:end-1);
    t = [ti t];
    ti = linspace(x(end), b, L+1);
    ti = ti(2:end-1);
    t = [t ti];
end
