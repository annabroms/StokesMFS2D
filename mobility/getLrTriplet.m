function Lr = getLrTriplet(B1,B2,B3,K1,K2,K3)
%GETLRTRIPLET Construct the triplet rigid-motion closure matrix.
%
% Input ordering:
%   [x_1; x_2; x_3; y_1; y_2; y_3]

if nargin < 6
    error('getLrTriplet requires B1, B2, B3, K1, K2, and K3.');
end

L1 = B1*K1';
L2 = B2*K2';
L3 = B3*K3';

[xx1,xy1,yx1,yy1] = split_block(L1,size(B1,1)/2,size(K1,1)/2);
[xx2,xy2,yx2,yy2] = split_block(L2,size(B2,1)/2,size(K2,1)/2);
[xx3,xy3,yx3,yy3] = split_block(L3,size(B3,1)/2,size(K3,1)/2);

Lxx = blkdiag(xx1,xx2,xx3);
Lxy = blkdiag(xy1,xy2,xy3);
Lyx = blkdiag(yx1,yx2,yx3);
Lyy = blkdiag(yy1,yy2,yy3);

Lr = [Lxx, Lxy; Lyx, Lyy];

end

function [xx,xy,yx,yy] = split_block(L,nrow,ncol)
xx = L(1:nrow,1:ncol);
xy = L(1:nrow,ncol+1:end);
yx = L(nrow+1:end,1:ncol);
yy = L(nrow+1:end,ncol+1:end);
end
