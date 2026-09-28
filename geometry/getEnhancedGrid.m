function [cent_clust_cells, acc_cells, coll_clust_cells, clust_pairs, coll_pairs, pairs] = getEnhancedGrid(q, opt)
%GETENHANCEDGRID Clustered ellipse-segment nodes for close particle pairs.
%   [cent_clust_cells, acc_cells, coll_clust_cells, clust_pairs, coll_pairs, pairs] = getEnhancedGrid(q, opt)
%
%   Inputs:
%     q   - complex centers (P x 1)
%     opt - struct with fields:
%           rad    : radii (P x 1), default ones
%           Nclust  : number of clustered nodes per close pair
%           beta    : ellipse tip parameter (0<beta<1)
%           r_proxy : proxy radius for filtering |z|>r_proxy
%           delta_pair : proximity threshold for pairs
%           ellipse_constant : if true, build the ellipse segment (and its
%                    r_proxy node filter) for every close pair using the
%                    fixed gap opt.smallest_delta instead of the pair's
%                    actual gap. Node placement still follows the true
%                    centers, so only the discretisation (which nodes pass
%                    the r_proxy filter) is frozen, removing the jumps that
%                    occur when nodes cross that filter as the true gap
%                    changes. Requires opt.smallest_delta.
%           smallest_delta : smallest gap that will ever be requested for a
%                    close pair; only used when opt.ellipse_constant=true.
%           visualise_grid : logical, optional plotting
%                    flag for node visualisation
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

if isfield(opt,'rad') && ~isempty(opt.rad)
    rad = opt.rad;
    if isscalar(rad)
        rad = repmat(rad,P,1);
    end
    if numel(rad) ~= P
        error('getEnhancedGrid:badRadii','opt.rad must be scalar or length P.');
    end
else
    rad = ones(P,1);
end


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
    delta_pair = 0.2;
end

ellipse_constant = isfield(opt,'ellipse_constant') && logical(opt.ellipse_constant);
if ellipse_constant
    if ~isfield(opt,'smallest_delta') || isempty(opt.smallest_delta)
        error('getEnhancedGrid:MissingSmallestDelta', ...
            'opt.ellipse_constant=true requires opt.smallest_delta.');
    end
    smallest_delta = opt.smallest_delta;
end

if isfield(opt,'visualise_grid')
    visualise_grid = logical(opt.visualise_grid);
else
    visualise_grid = false;
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
        ri = rad(i); rj = rad(j);
        D = abs(cj-ci);
        gap = D - (ri + rj);
        if gap < delta_pair
            pairs = [pairs; i j];

            if ellipse_constant
                if (smallest_delta-gap)>1e-9;
                    error('getEnhancedGrid:GapBelowSmallestDelta', ...
                        ['Pair (%d,%d) has gap=%.3g, smaller than ', ...
                         'opt.smallest_delta=%.3g.'],i,j,gap,smallest_delta);
                end
                gap_ellipse = smallest_delta;
            else
                gap_ellipse = gap;
            end

            [cent_i, cent_j, coll_i, coll_j, zacc_i, zacc_j] = ...
                pair_clusters_ellipse(ci, cj, ri, rj, Nclust, gap_ellipse, r_proxy, beta);

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

if visualise_grid
    showEnhancedGridPoints(q,rad,r_proxy,cent_clust_cells,acc_cells, ...
        coll_clust_cells,clust_pairs,coll_pairs,pairs);
end
end

function getEnhancedGrid_selftest()
% Self-test for three particles with uniform proxy circles
P = 3;
q = [0; 2.1; 1.05 + 1.8i];
opt.rad = ones(P,1);
opt.Nclust = 100;
opt.beta = 0.3;
opt.r_proxy = 0.7;
opt.delta_pair = 0.2;
opt.visualise_grid = true;

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
plot(proxy, 'r--'); hold on;
%plot(q, 'k*', 'MarkerSize', 8);
plot(acc_pts, 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'y');
plot(cent_clust, 'r--');
plot(coll_clust, 'b.');
axis equal; grid on;
legend('proxy boundaries','centers','accumulation points','clustered centers','clustered collocation','Location','best');
end

function showEnhancedGridPoints(q,rad,r_proxy,cent_clust_cells,acc_cells, ...
    coll_clust_cells,clust_pairs,coll_pairs,pairs)
%SHOWENHANCEDGRIDPOINTS Plot all node families used in getEnhancedGrid.

draw_fine_grid = 0; % visualise fine discretization. 

P = numel(q);
t = linspace(0,2*pi,240).';

figure('Name','getEnhancedGrid points');
hold on;

for k = 1:P
    body = q(k) + rad(k)*(cos(t)+1i*sin(t));
    proxy = q(k) + r_proxy*(cos(t)+1i*sin(t));
    if k == 1
        fill(real(body), imag(body), ...
        [1 0.5 0], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.2,'DisplayName','Active particle');


        plot(real(body),imag(body),'b.','MarkerSize',10,'DisplayName','Coarse collocation nodes');
        plot(real(proxy),imag(proxy),'r.','DisplayName','Fine source points'); % $\mathbf{\mathcal{Y}}^{(1\text{-}2)}$');
    else
       % plot(real(body),imag(body),'k-','LineWidth',1.0,'HandleVisibility','off');
        plot(real(proxy),imag(proxy),'r.','HandleVisibility','off');
        
    end
    if draw_fine_grid
        plot(real(body),imag(body),'k-','HandleVisibility','off');
    end
end

%plot(real(q),imag(q),'kp','MarkerSize',9,'MarkerFaceColor','k','DisplayName','particle centers');
% for k = 1:P
%     text(real(q(k)),imag(q(k)),sprintf('  q_%d',k), ...
%         'Color','k','FontSize',10,'Interpreter','none');
% end

for k = 1:P
    rk = cent_clust_cells{k};
    if ~isempty(rk)
        ind = abs(rk-q(k))<r_proxy;
        if k == 1
            plot(real(rk(ind)),imag(rk(ind)),'.','Color',[1 0 1 0.1],'MarkerSize',2,'DisplayName','Discarded ellipse points');
        else
            plot(real(rk),imag(rk),'.','Color',[1 0 1 0.1],'MarkerSize',2,'HandleVisibility','off');
        end
        ind = abs(rk-q(k))>r_proxy;
        plot(real(rk(ind)),imag(rk(ind)),'r.','HandleVisibility','off');
        if draw_fine_grid
            plot(real([q(k), acc_cells{k}]),imag([q(k), acc_cells{k}]),'k--','HandleVisibility','off');
            plot(real(q(k)),imag(q(k)),'k.','MarkerSize',10,'HandleVisibility','off');
        end
    end
end

for k = 1:P
    ck = coll_clust_cells{k};
    if ~isempty(ck)
         if k == 1
             plot(real(ck),imag(ck),'m.','MarkerSize',8,'DisplayName','Fine collocation nodes');
         else
             plot(real(ck),imag(ck),'m.','MarkerSize',8,'HandleVisibility','off');
         end
         %   plot(real(ck),imag(ck),'b.','MarkerSize',8,'HandleVisibility','off');
 %       end
    end
end


xlabel('x');
ylabel('y');
%title('Enhanced-grid node families and pair labels','Interpreter','none');
%legend('Interpreter','latex','FontSize',16);

if draw_fine_grid

for row = 1:size(pairs,1)
    i = pairs(row,1);
    j = pairs(row,2);
    plot(real([q(i); q(j)]),imag([q(i); q(j)]),'Color',[0.25 0.25 0.25], ...
        'LineStyle',':','LineWidth',1.0,'HandleVisibility','off');
end

for k = 1:P
    ak = acc_cells{k};
    if ~isempty(ak)
        if k == 1
            plot(real(ak),imag(ak),'b.','MarkerSize',10,'MarkerFaceColor','y', ...
                'DisplayName','Image accumulation points');
        else
            plot(real(ak),imag(ak),'b.','MarkerSize',10,'MarkerFaceColor','y', ...
                'HandleVisibility','off');
        end
    end
end

x1 = q(2);
x2 = acc_cells{2};
y1 = -0.2;
y2 = y1;
draw_arrow(x1,y1,x2,y2,0.2);

%%
x2 = q(2)-1;
y1 = y1-0.4;
y2 = y1; 
draw_arrow(x1,y1,x2,y2,0.2);

%%

x1 = q(2);
y1 = 0;
y2 = r_proxy*sin(pi/5);
x2 = q(2)+r_proxy*cos(pi/5);

draw_arrow(x1,y1,x2,y2,0.3);

%%
plot([q(2)-1,q(2)-1],[-0.6,0.5],'k-','HandleVisibility','off')
plot([q(1)+1,q(1)+1],[-0.6,0.5],'k-','HandleVisibility','off')
%%
% x2 = q(1)+1;
% x1 = q(1)+0.92;
% y1 = 0.46;
% y2 = 0.46;
% plot([x1 x2], [y1 y2], 'k-', 'LineWidth', 1.5,'HandleVisibility','off')
% 
% % Direction vector
% dx = x2 - x1;
% dy = y2 - y1;
% 
% % Arrowheads at both ends
% quiver(x1, y1,  dx,  dy, 0, 'k', 'LineWidth', 1.3, 'MaxHeadSize', 70,'HandleVisibility','off')

%%
x1 = q(2)-1;
x2 = q(2)-0.92;
y1 = 0.46;
y2 = 0.46;

% plot([x1 x2], [y1 y2], 'k-', 'LineWidth', 1.5,'HandleVisibility','off')
% 
% % Direction vector
% dx = x2 - x1;
% dy = y2 - y1;
% 
% % Arrowheads at both ends
% quiver(x2, y2, -dx, -dy, 0, 'k', 'LineWidth', 1.3, 'MaxHeadSize', 70,'HandleVisibility','off')



%plot([q(1)+1,q(2)-1],[0.55,0.55],'k-','LineWidth', 1.5,'HandleVisibility','off')

end

legend show
set(legend,'interpreter','latex','FontSize',20,'box','off')

axis equal;
axis off
end

function draw_arrow(x1,y1,x2,y2,aSize)

plot([x1 x2], [y1 y2], 'k-', 'LineWidth', 1.5,'HandleVisibility','off')

% Direction vector
dx = x2 - x1;
dy = y2 - y1;

% Arrowheads at both ends
quiver(x1, y1,  dx,  dy, 0, 'k', 'LineWidth', 1.3, 'MaxHeadSize', aSize,'HandleVisibility','off')
quiver(x2, y2, -dx, -dy, 0, 'k', 'LineWidth', 1.3, 'MaxHeadSize', aSize,'HandleVisibility','off')

end
