function Kft_triplet = getKftTriplet(K1,K2,K3)
%GETKFTTRIPLET Map triplet-ordered source strengths to per-body RBM moments.
%
% Input ordering:
%   [f_x1; f_x2; f_x3; f_y1; f_y2; f_y3]
%
% Output ordering:
%   [F_x1; F_y1; T_1; F_x2; F_y2; T_2; F_x3; F_y3; T_3]

if nargin < 3
    error('getKftTriplet requires K1, K2, and K3.');
end

[K1x,K1y] = split_k(K1);
[K2x,K2y] = split_k(K2);
[K3x,K3y] = split_k(K3);

n1 = size(K1x,2);
n2 = size(K2x,2);
n3 = size(K3x,2);

Kft_triplet = [ ...
    K1x, zeros(3,n2), zeros(3,n3), K1y, zeros(3,n2), zeros(3,n3); ...
    zeros(3,n1), K2x, zeros(3,n3), zeros(3,n1), K2y, zeros(3,n3); ...
    zeros(3,n1), zeros(3,n2), K3x, zeros(3,n1), zeros(3,n2), K3y];

end

function [Kx,Ky] = split_k(K)
Kt = K';
n = size(K,1)/2;
Kx = Kt(:,1:n);
Ky = Kt(:,n+1:end);
end
