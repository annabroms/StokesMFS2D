clear;
close all;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Random MC Close Neighbors (Apr 14, 2026) ===\n');

phi = 0.55;
domain = 'boxed';
rad = 1;
min_gap = 1e-3;
delta_pair = 0.2;
n_sweeps = 20;

P_values = 20:20:500;
P_values = 50:50:1000;
nruns = 4;
seed0 = 100;

nP = numel(P_values);
close_pair_counts = nan(nP,nruns);
close_neighbors_per_particle = nan(nP,nruns);
min_gaps = nan(nP,nruns);

fprintf('phi = %.3f, domain = %s, delta_pair = %.3f, nruns = %d\n', ...
    phi,domain,delta_pair,nruns);

for ip = 1:nP
    P = P_values(ip);
    fprintf('\nP = %d\n',P);

    for irun = 1:nruns
        geom_opt = struct();
        geom_opt.domain = domain;
        geom_opt.phi = phi;
        geom_opt.rad = rad;
        geom_opt.min_gap = min_gap;
        geom_opt.n_sweeps = n_sweeps;
        geom_opt.rng_seed = seed0 + 1000*ip + irun;
        geom_opt.visualise = false;

        [q, meta] = random_discs_mc(P, geom_opt);
        n_close = count_close_pairs(q, delta_pair, rad);

        close_pair_counts(ip,irun) = n_close;
        close_neighbors_per_particle(ip,irun) = 2*n_close / P;
        min_gaps(ip,irun) = meta.min_surface_gap;

        fprintf('  run %d/%d: close pairs = %d, close neighbors/particle = %.3f, min gap = %.3e\n', ...
            irun,nruns,n_close,close_neighbors_per_particle(ip,irun),min_gaps(ip,irun));
    end
end

pair_mean = mean(close_pair_counts,2);
pair_min = min(close_pair_counts,[],2);
pair_max = max(close_pair_counts,[],2);

neigh_mean = mean(close_neighbors_per_particle,2);
neigh_min = min(close_neighbors_per_particle,[],2);
neigh_max = max(close_neighbors_per_particle,[],2);

fprintf('\nSummary:\n');
fprintf('  %-8s %-12s %-12s %-12s %-12s\n', ...
    'P','pair_min','pair_max','pair_mean','neigh_mean');
for ip = 1:nP
    fprintf('  %-8d %-12.0f %-12.0f %-12.3f %-12.3f\n', ...
        P_values(ip),pair_min(ip),pair_max(ip),pair_mean(ip),neigh_mean(ip));
end

figure('Name','apr14 random mc close neighbors','Color','w');
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
hold on;
for irun = 1:nruns
    plot(P_values,close_pair_counts(:,irun),'o','Color',[0.75 0.75 0.75], ...
        'MarkerFaceColor',[0.85 0.85 0.85], ...
        'DisplayName',ternary(irun == 1,'individual runs',''));
end
plot_band(P_values,pair_min,pair_max,[0.00 0.45 0.74]);
plot(P_values,pair_mean,'-o','Color',[0.00 0.45 0.74], ...
    'LineWidth',1.5, ...
    'MarkerFaceColor',[0.00 0.45 0.74], ...
    'DisplayName','mean close pairs');
grid on;
xlabel('P');
ylabel('Number of close pairs');
title(sprintf('Close pairs with gap < %.3f at \\phi = %.3f',delta_pair,phi), ...
    'Interpreter','none');
legend('Location','northwest');

nexttile;
hold on;
plot_band(P_values,neigh_min,neigh_max,[0.85 0.33 0.10]);
plot(P_values,neigh_mean,'-o','Color',[0.85 0.33 0.10], ...
    'LineWidth',1.5, ...
    'MarkerFaceColor',[0.85 0.33 0.10], ...
    'DisplayName','mean close neighbors / particle');
grid on;
xlabel('P');
ylabel('Close neighbors per particle');
title('Average close-neighbor count per particle','Interpreter','none');
legend('Location','northwest');

title(t,'Random MC close-neighbor statistics','Interpreter','none');

function plot_band(x,ymin,ymax,color)
xx = [x(:); flipud(x(:))];
yy = [ymin(:); flipud(ymax(:))];
patch(xx,yy,color, ...
    'FaceAlpha',0.15, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');
end

function out = ternary(cond,a,b)
if cond
    out = a;
else
    out = b;
end
end
