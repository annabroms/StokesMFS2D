%PARTICLE_LINE_ELAST Elastance solve for a line of particles.
%
% This demo increases the number of particles P in a straight line with a
% fixed particle-particle gap delta, solves the elastance problem with
% the peanut-enhanced solver, and records the maximum relative boundary
% residual and GMRES iteration count. If capacitance = true, it also solves
% the matching capacitance problem on the same geometries and reports those
% results in separate figures.

clear;
close all;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

rad = 2;
delta_vec = [1e-1 1e-2 1e-3];
delta_labels = arrayfun(@(d) sprintf('$\\delta = %.3g$', d), ...
    delta_vec, 'UniformOutput', false);
P_values = 50:20:400;
font_size = 14;  % Controls axes, labels, titles, and legends.

N_peanut = 400;
%delta_pair = 0.2;
delta_pair_cap = 2*rad+0.6; % works much better for capacitance! The pairs are coupled. WHY?
delta_pair_elast = 0.2;
gmres_tol = 1e-7;
capacitance = true;

max_res_values = zeros(length(delta_vec),numel(P_values));
it_values = zeros(length(delta_vec),numel(P_values));
if capacitance
    cap_max_res_values = zeros(length(delta_vec),numel(P_values));
    cap_it_values = zeros(length(delta_vec),numel(P_values));
end

fprintf('=== particle_line_elast elastance demo ===\n');

N_c = 60;
N_f = 60; 
opt = getLaplace2Dparams(1,rad,N_c,N_f);
opt.N_peanut = N_peanut;
opt.delta_pair = delta_pair_elast;
opt.gmres_tol = gmres_tol;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.debug = 0;
opt.get_bndry_field = true;
opt.reuse_pair_basis_by_sep = true;
opt_cap = opt;
opt_cap.delta_pair = delta_pair_cap;

for i = 1:length(delta_vec)
    delta = delta_vec(i);

    for k = 1:length(P_values)
        P = P_values(k);
        q = buildParticleLine(P,rad,delta);
        v_body = buildAlternatingVoltages(q,rad);

        opt.P = P;
      
        fprintf('Solving delta = %.3g, P=%u ...\n',delta,P);
        [~,sol] = solve_elast_peanut(q,v_body,opt);

        max_res_values(i,k) = sol.maxres;
        it_values(i,k) = sol.it;
        fprintf('  maxres = %.3e, it = %d\n',max_res_values(i,k),it_values(i,k));

        if capacitance
            fprintf('Solving capacitance delta = %.3g, P=%u ...\n',delta,P);
            [~,cap_sol] = solve_cap_peanut(q,v_body,opt_cap);

            cap_max_res_values(i,k) = cap_sol.maxres;
            cap_it_values(i,k) = cap_sol.it;
            fprintf('  maxres = %.3e, it = %d\n',cap_max_res_values(i,k),cap_it_values(i,k));
        end
    end
end

%%
figure('Name','Elastance residual vs P');
for i = 1:length(delta_vec)
    semilogy(P_values,max_res_values(i,:),'-o', ...
        'LineWidth',1.5,'MarkerSize',6, ...
        'DisplayName',delta_labels{i});
    hold on
end
styleLinePlot(gca,font_size,'Number of particles $P$', ...
    'Maximum boundary residual', ...
    'Elastance peanut solve residual vs $P$');

figure('Name','Elastance iterations vs P');
hold on;
for i = 1:length(delta_vec)
    plot(P_values,it_values(i,:),'-o', ...
        'LineWidth',1.5,'MarkerSize',6, ...
        'DisplayName',delta_labels{i});
end
styleLinePlot(gca,font_size,'Number of particles $P$', ...
    'GMRES iterations', ...
    'Elastance peanut solve iterations vs $P$');

if capacitance
    figure('Name','Capacitance residual vs P');
    hold on;
    for i = 1:length(delta_vec)
        semilogy(P_values,cap_max_res_values(i,:),'-o', ...
        'LineWidth',1.5,'MarkerSize',6, ...
        'DisplayName',delta_labels{i});
    end
    styleLinePlot(gca,font_size,'Number of particles $P$', ...
        'Maximum boundary residual', ...
        'Capacitance peanut solve residual vs $P$');

    figure('Name','Capacitance iterations vs P');
    hold on;
    for i = 1:length(delta_vec)
        plot(P_values,cap_it_values(i,:),'-o', ...
        'LineWidth',1.5,'MarkerSize',6, ...
        'DisplayName',delta_labels{i});
    end
    styleLinePlot(gca,font_size,'Number of particles $P$', ...
        'GMRES iterations', ...
        'Capacitance peanut solve iterations vs $P$');
end

function q = buildParticleLine(P,rad,delta)
    spacing = 2*rad + delta*rad;
    offset = (0:P-1)' - (P-1)/2;
    q = spacing*offset;
    q = complex(q,zeros(size(q)));
end
