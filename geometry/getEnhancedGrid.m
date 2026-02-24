function [cent_clust_cells, acc_cells, coll_clust_cells, clust_pairs, coll_pairs, pairs] = getEnhancedGrid(q, opt)
%GETENHANCEDGRID Clustered ellipse-segment nodes for close particle pairs.
%   [cent_clust_cells, acc_cells, coll_clust_cells, clust_pairs, coll_pairs, pairs] = getEnhancedGrid(q, opt)
%
%   Inputs:
%     q   - complex centers (P x 1)
%     opt - struct with fields:
%           rads    : radii (P x 1), default ones
%           Nclust  : number of clustered nodes per close pair
%           beta    : ellipse tip parameter (0<beta<1)
%           r_proxy : proxy radius for filtering |z|>r_proxy
%           delta_pair : proximity threshold for pairs
%
%   Outputs:
%     cent_clust_cells : cell(P,1) of clustered centers for each particle
%     acc_cells        : cell(P,1) of accumulation points per particle
%     coll_clust_cells : cell(P,1) of clustered collocation nodes per particle
%     clust_pairs      : cell(P,P) of clustered centers per close pair
%     coll_pairs       : cell(P,P) of clustered collocation nodes per close pair
%     pairs            : N x 2 array of close pairs [i j] with i<j
%
%   Self-test (no inputs): visualises discretization and proxy circles for P=3.

if nargin == 0
    getEnhancedGrid_selftest();
    return;
end

q = q(:);
P = numel(q);

rads = ones(P,1);
Nclust = opt.Nclust;
beta = opt.beta;
if isfield(opt,'Rp_f')
    r_proxy = opt.Rp_f;
elseif isfield(opt,'r_proxy')
    r_proxy = opt.r_proxy;
else
    error('getEnhancedGrid:missingProxyRadius','Need opt.Rp_f or opt.r_proxy.');
end
if isfield(opt, 'delta_pair')
    delta_pair = opt.delta_pair;
else
    delta_pair = 5e-2;
end

cent_clust_cells = cell(P,1);
coll_clust_cells = cell(P,1);
acc_cells = cell(P,1);
clust_pairs = cell(P,P);
coll_pairs = cell(P,P);
pairs = [];

for i = 1:P-1
    for j = i+1:P
        ci = q(i); cj = q(j);
        ri = rads(i); rj = rads(j);
        D = abs(cj-ci);
        gap = D - (ri + rj);
        if gap < delta_pair
            pairs = [pairs; i j];

            [cent_i, cent_j, coll_i, coll_j, zacc_i, zacc_j] = ...
                pair_clusters_ellipse(ci, cj, ri, rj, Nclust, gap, r_proxy, beta);

            % Add enhancement nodes only if the accumulation point lies
            % outside the proxy radius for that particle.
            add_i = abs(zacc_i - ci) > r_proxy;
            add_j = abs(zacc_j - cj) > r_proxy;

            if add_i
                cent_clust_cells{i} = [cent_clust_cells{i}; cent_i];
                coll_clust_cells{i} = [coll_clust_cells{i}; coll_i];
                acc_cells{i} = [acc_cells{i}; zacc_i];
                clust_pairs{i,j} = cent_i;
                coll_pairs{i,j} = coll_i;
            end

            if add_j
                cent_clust_cells{j} = [cent_clust_cells{j}; cent_j];
                coll_clust_cells{j} = [coll_clust_cells{j}; coll_j];
                acc_cells{j} = [acc_cells{j}; zacc_j];
                clust_pairs{j,i} = cent_j;
                coll_pairs{j,i} = coll_j;
            end
        end
    end
end
end

function getEnhancedGrid_selftest()
% Self-test for three particles with uniform proxy circles
P = 3;
q = [0; 2.1; 1.05 + 1.8i];
opt.rads = ones(P,1);
opt.Nclust = 100;
opt.beta = 0.5;
opt.r_proxy = 0.7;
opt.delta_pair = 0.5;

[cent_clust_cells, acc_cells, coll_clust_cells, clust_pairs, coll_pairs, pairs] = getEnhancedGrid(q, opt);
cent_clust = vertcat(cent_clust_cells{:});
acc_pts = vertcat(acc_cells{:});
coll_clust = vertcat(coll_clust_cells{:});

% Proxy boundaries
N = 200;
t = linspace(0,2*pi,N).';
proxy = [];
for k = 1:P
    proxy = [proxy; q(k) + opt.r_proxy * (cos(t) + 1i*sin(t))];
end

figure('Name','getEnhancedGrid selftest');
plot(proxy, 'k--'); hold on;
plot(q, 'k*', 'MarkerSize', 8);
plot(acc_pts, 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'y');
plot(cent_clust, 'ro');
plot(coll_clust, 'b.');
axis equal; grid on;
legend('proxy boundaries','centers','accumulation points','clustered centers','clustered collocation','Location','best');
end
