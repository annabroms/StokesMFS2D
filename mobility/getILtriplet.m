function ILtriplet = getILtriplet(L)
%GETILTRIPLET Construct the triplet analogue of getILpair.
%
% Input ordering:
%   [x_1; x_2; x_3; y_1; y_2; y_3]

if nargin < 1
    error('getILtriplet requires a single-body projection matrix L.');
end

n = size(L,1)/2;
Lxx = L(1:n,1:n);
Lxy = L(1:n,n+1:end);
Lyx = L(n+1:end,1:n);
Lyy = L(n+1:end,n+1:end);

Axx = blkdiag(eye(n)-Lxx, eye(n)-Lxx, eye(n)-Lxx);
Axy = blkdiag(-Lxy, -Lxy, -Lxy);
Ayx = blkdiag(-Lyx, -Lyx, -Lyx);
Ayy = blkdiag(eye(n)-Lyy, eye(n)-Lyy, eye(n)-Lyy);

ILtriplet = [Axx, Axy; Ayx, Ayy];

end
