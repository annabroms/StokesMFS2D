function [cent_i, cent_j, coll_i, coll_j, zacc_i, zacc_j] = pair_clusters_ellipse(ci, cj, ri, rj, Nclust, gap, r_proxy, beta)
%PAIR_CLUSTERS_ELLIPSE Clustered source and collocation nodes for a close pair.
%   Uses analytic accumulation points and ellipse-segment nodes.

D = abs(cj-ci);
if D == 0
    cent_i = []; cent_j = []; coll_i = []; coll_j = [];
    zacc_i = []; zacc_j = [];
    return;
end

% Analytic accumulation points (circles)
[acc1, acc2] = accumulation_points(ri, rj, gap, [real(ci); imag(ci)], [real(cj); imag(cj)]);
zacc_i = acc1(1) + 1i*acc1(2);
zacc_j = acc2(1) + 1i*acc2(2);

v = (cj-ci)/D;

% Ellipse segments (shifted Chebyshev)
acc1n = abs(zacc_i - ci) / ri;
acc2n = abs(zacc_j - cj) / rj;
[ell1, ~] = ellipse_cheb_segment(acc1n, Nclust, r_proxy/ri, beta);
[ell2, ~] = ellipse_cheb_segment(acc2n, Nclust, r_proxy/rj, beta);

cent_i = ci + (ri*ell1) * v;
cent_j = cj + (rj*ell2) * (-v);

% Extra collocation clustered near closest approach
Mclust = max(50, 5*Nclust);
t = linspace(0,2*pi,Mclust).';
tclust = pi + mobius_map(t, 0.9); %0.7

ang_i = angle(v);
ang_j = angle(-v);
coll_i = ci + ri*(cos(ang_i + tclust) + 1i*sin(ang_i + tclust));
coll_j = cj + rj*(cos(ang_j + tclust) + 1i*sin(ang_j + tclust));
end

function theta = mobius_map(t, r)
%MOBIUS_MAP  Cluster angles near theta=0 via Möbius map.
    z = exp(1i*t);
    w = (z - r) ./ (1 - r*z);
    theta = mod(angle(w), 2*pi);
end
