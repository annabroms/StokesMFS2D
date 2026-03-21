%PARTICLE_LINE_MOB Mobility solve for a line of particles.
%
% This demo increases the number of particles P in a straight line with a
% fixed particle-particle gap delta, solves the mobility problem with the
% peanut-enhanced solver, and records the relative boundary residual and
% GMRES iteration count. If resistance = true, it also solves the matching
% resistance problem on the same geometries and reports those results in
% separate figures.

clear;
close all;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

rad = 1;
delta_vec = [1e-1 1e-2 1e-3];
delta_labels = arrayfun(@(d) sprintf('$\\delta = %.3g$', d), ...
    delta_vec, 'UniformOutput', false);
P_values = 10:10:100;

N_c = 100;
N_peanut = 400;
delta_pair = 0.2;
gmres_tol = 1e-8;
load_mode = 'alternating';  % 'smooth' keeps the original trigonometric loads.
resistance = true;

fprintf('=== particle_line mobility demo ===\n');

opt = get2Dparams(1,N_c);
opt.N_peanut = N_peanut;
opt.rad = rad;
opt.delta_pair = delta_pair;
opt.gmres_tol = gmres_tol;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.debug = 0;
opt.get_bndry_field = true;
opt.reuse_pair_basis_by_sep = true;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt_res = opt;

rel_res_values = zeros(numel(delta_vec),numel(P_values));
it_values = zeros(numel(delta_vec),numel(P_values));
if resistance
    res_rel_res_values = zeros(numel(delta_vec),numel(P_values));
    res_it_values = zeros(numel(delta_vec),numel(P_values));
end

for i = 1:length(delta_vec)
    delta = delta_vec(i);
    for k = 1:length(P_values)
        P = P_values(k);
        opt.P = P;
    
        % Particles along a line
        q = (2*rad + delta)*(0:P-1)';
        q = complex(q,zeros(size(q)));
        [F,T] = buildLoads(q,load_mode);
    
        fprintf('Solving delta = %.3g, P=%u ...\n',delta,P);
        [~,sol] = solve_mob_peanut_enhanced(q,F,T,opt);
    
        rel_res_values(i,k) = sol.rel_res;
        it_values(i,k) = sol.it;
        fprintf('  rel_res = %.3e, it = %d\n',rel_res_values(i,k),it_values(i,k));

        if resistance
            [U,W] = buildResistanceMotion(q,load_mode);
            fprintf('Solving resistance delta = %.3g, P=%u ...\n',delta,P);
            [~,res_sol] = solve_res_peanut_enhanced(q,U,W,opt_res);

            res_rel_res_values(i,k) = res_sol.rel_res;
            res_it_values(i,k) = res_sol.it;
            fprintf('  rel_res = %.3e, it = %d\n',res_rel_res_values(i,k),res_it_values(i,k));
        end
    end

end

figure('Name','Mobility residual vs P');
hold on;
for i = 1:length(delta_vec)
    semilogy(P_values,rel_res_values(i,:),'-o', ...
        'LineWidth',1.5,'MarkerSize',6, ...
        'DisplayName',delta_labels{i});
end
xlabel('Number of particles P');
ylabel('Relative boundary residual');
title('Mobility peanut solve residual vs P');
grid on;
axis tight;
legend('Location','best','Interpreter','latex');


figure('Name','Mobility iterations vs P');
hold on;
for i = 1:length(delta_vec)
    plot(P_values,it_values(i,:),'-o', ...
        'LineWidth',1.5,'MarkerSize',6, ...
        'DisplayName',delta_labels{i});
end
xlabel('Number of particles P');
ylabel('GMRES iterations');
title('Mobility peanut solve iterations vs P');
grid on;
axis tight;
legend('Location','best','Interpreter','latex');

if resistance
    figure('Name','Resistance residual vs P');
    hold on;
    for i = 1:length(delta_vec)
        semilogy(P_values,res_rel_res_values(i,:),'-o', ...
            'LineWidth',1.5,'MarkerSize',6, ...
            'DisplayName',delta_labels{i});
    end
    xlabel('Number of particles P');
    ylabel('Relative boundary residual');
    title('Resistance peanut solve residual vs P');
    grid on;
    axis tight;
    legend('Location','best','Interpreter','latex');

    figure('Name','Resistance iterations vs P');
    hold on;
    for i = 1:length(delta_vec)
        plot(P_values,res_it_values(i,:),'-o', ...
            'LineWidth',1.5,'MarkerSize',6, ...
            'DisplayName',delta_labels{i});
    end
    xlabel('Number of particles P');
    ylabel('GMRES iterations');
    title('Resistance peanut solve iterations vs P');
    grid on;
    axis tight;
    legend('Location','best','Interpreter','latex');
end


function [F,T] = buildLoads(q,load_mode)
x = real(q);
n = numel(x);

switch load_mode
    case 'smooth'
        F = [cos(0.3*x), sin(0.5*x)];
        T = 0.25*cos(0.7*x);
    case 'alternating'
        s = ones(n,1);
        s(2:2:end) = -1;
        F = [s, zeros(n,1)];
        T = s;
    otherwise
        error('Unknown load_mode: %s',load_mode);
end
end

function [U,W] = buildResistanceMotion(q,load_mode)
x = real(q);
n = numel(x);

switch load_mode
    case 'smooth'
        U = [cos(0.3*x), sin(0.5*x)];
        W = 0.25*cos(0.7*x);
    case 'alternating'
        s = ones(n,1);
        s(2:2:end) = -1;
        U = [s, zeros(n,1)];
        W = s;
    otherwise
        error('Unknown load_mode: %s',load_mode);
end
end
