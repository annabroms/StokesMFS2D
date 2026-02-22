function [z, theta, alpha, c, s] = ellipse_cheb_segment(acc, M, r_proxy, beta)
%ELLIPSE_CHEB_SEGMENT Chebyshev nodes on ellipse with focus at acc, filtered by |z|>r_proxy.
%   [z, theta, alpha, c, s] = ellipse_cheb_segment(acc, N, r_proxy)
%   assumes 0 < acc < 1 (real). Nodes are generated from
%       z = c + s*cos(theta + i*alpha)
%   with foci at c ± s = {0, acc}, and alpha chosen so that the rightmost tip
%   is at acc + beta*(1-acc). Points with |z| <= r_proxy are discarded.

if acc <= 0 || acc >= 1
    error('acc must satisfy 0 < acc < 1.');
end
if M < 1
    z = [];
    theta = [];
    alpha = [];
    c = [];
    s = [];
    %error('N must be positive.');
end
if nargin < 3
    r_proxy = 0;
end
if nargin < 4
    beta = 0.5;
end
if beta <= 0 || beta >= 1
    error('beta must satisfy 0 < beta < 1.');
end

c = acc/2;
s = acc/2;
tip = acc + beta*(1-acc);
alpha = acosh((tip - c) / s);

j = (0:2*M-1).';
theta = (j + 0.5) * pi / M;
z_all = c + s * cos(theta + 1i*alpha);
keep = abs(z_all) > r_proxy;
z = z_all(keep);
theta = theta(keep);   

end
