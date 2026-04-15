function Lf = getLfTriplet(K1,K2,K3)
%GETLFTRIPLET Construct the triplet fine-source rigid projection matrix.
%
% Input ordering:
%   [x_1; x_2; x_3; y_1; y_2; y_3]

if nargin < 3
    error('getLfTriplet requires K1, K2, and K3.');
end

L1 = K1*((K1'*K1)\K1');
L2 = K2*((K2'*K2)\K2');
L3 = K3*((K3'*K3)\K3');

[xx1,xy1,yx1,yy1] = split_block(L1);
[xx2,xy2,yx2,yy2] = split_block(L2);
[xx3,xy3,yx3,yy3] = split_block(L3);

Lxx = blkdiag(xx1,xx2,xx3);
Lxy = blkdiag(xy1,xy2,xy3);
Lyx = blkdiag(yx1,yx2,yx3);
Lyy = blkdiag(yy1,yy2,yy3);

Lf = [Lxx, Lxy; Lyx, Lyy];

end

function [xx,xy,yx,yy] = split_block(L)
n = size(L,1)/2;
xx = L(1:n,1:n);
xy = L(1:n,n+1:end);
yx = L(n+1:end,1:n);
yy = L(n+1:end,n+1:end);
end
