close all;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Random Discs Peanut Mobility (Apr 14, 2026) ===\n');

% Geometry
geom.P = 200;
geom.rad = 1;
geom.domain = 'boxed';
geom.phi = 0.65;
geom.min_gap = 1e-3;

% Monte Carlo generation
mc.n_sweeps = 30;
mc.rng_seed = 1;
mc.show_generation_plot = false;

% Applied loads
loads.rng_seed = 11;

% Solver
solver.N_c = 60;
solver.N_f = 150;
solver.N_peanut = 400;
solver.delta_pair = 0.2;
solver.gmres_tol = 1e-8;
solver.maxit = 1000;

% Reporting and visualisation
report.count_close_pairs = (geom.P <= 3000) && strcmp(geom.domain,'boxed');
plots.make_figure = true;
plots.font_size = 16;
plots.full_edge_color = [0.10 0.10 0.10];
plots.line_width_full = 0.15;
plots.force_colormap = turbo(256);
plots.speed_colormap = parula(256);
plots.residual_colormap = hot(256);
plots.force_clim = [];
plots.speed_clim = [];
plots.residual_clim = [];

if geom.P <= 300
    plots.disk_points_full = 56;
elseif geom.P <= 1500
    plots.disk_points_full = 36;
elseif geom.P <= 5000
    plots.disk_points_full = 24;
else
    plots.disk_points_full = 18;
end

geom_opt = struct();
geom_opt.domain = geom.domain;
geom_opt.phi = geom.phi;
geom_opt.rad = geom.rad;
geom_opt.min_gap = geom.min_gap;
geom_opt.n_sweeps = mc.n_sweeps;
geom_opt.rng_seed = mc.rng_seed;
geom_opt.visualise = mc.show_generation_plot;

[q, geom_meta] = random_discs_mc(geom.P, geom_opt);

n_close = nan;
if report.count_close_pairs
    n_close = count_close_pairs(q, solver.delta_pair, geom.rad);
elseif strcmp(geom.domain,'periodic')
    fprintf('Close-pair count skipped: count_close_pairs is not periodic-image aware.\n');
else
    fprintf('Close-pair count skipped.\n');
end

opt = get2Dparams(geom.P, solver.N_c, solver.N_f);
opt.rad = geom.rad;
opt.delta_pair = solver.delta_pair;
opt.N_peanut = solver.N_peanut;
opt.gmres_tol = solver.gmres_tol;
opt.maxit = solver.maxit;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.debug = 0;
opt.gmres_verbose = 0;
opt.surface_error_mode = 'rel';
opt.reuse_pair_basis_by_sep = false;
opt.show_counter = 1;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.get_bndry_field = 1;

rng(loads.rng_seed,'twister');
F = randn(geom.P,2);
F = F - mean(F,1);
T = randn(geom.P,1);
T = ones(geom.P,1);
F = zeros(geom.P,2);

tic;
[UW, sol] = solve_mob_peanut_enhanced(q, F, T, opt);
solve_time = toc;

[U, W] = unpack_UW(UW);
body_forcing = sqrt(sum(F.^2,2) + (T/geom.rad).^2);
body_speed = sqrt(sum(U.^2,2)) + geom.rad*abs(W);
body_residual = sol.body_rel_res_max;
gmres_residual = final_gmres_residual(sol);

fprintf('Geometry:\n');
fprintf('  domain=%s, P=%d, phi_target=%.6f, phi=%.6f, L=%.6f\n', ...
    geom_meta.domain, geom.P, geom_meta.phi_target, geom_meta.phi, geom_meta.L);
fprintf('  min allowed gap=%.3e, measured min gap=%.3e\n', ...
    geom_meta.min_gap, geom_meta.min_surface_gap);
if report.count_close_pairs
    fprintf('  close pairs below delta_pair=%.3f: %d\n', solver.delta_pair, n_close);
end

fprintf('Solver:\n');
fprintf('  N_c=%d, N_f=%d, N_peanut=%d, gmres_tol=%.1e\n', ...
    opt.N_c, opt.N_f, opt.N_peanut, opt.gmres_tol);
fprintf('  iterations=%d, unknowns=%d, gmres_res=%.3e, surf_rel=%.3e, time=%.2fs\n', ...
    sol.it, sol.gmres_unknowns, gmres_residual, sol.rel_res, solve_time);

if isempty(body_residual)
    error('apr14_random_discs_peanut_mobility:MissingResiduals', ...
        'solve_mob_peanut_enhanced did not return sol.body_rel_res_max.');
end

max_body_residual = max(body_residual);
fprintf('  max(body_rel_res_max)=%.3e\n', max_body_residual);
fprintf('  |max(body_rel_res_max) - sol.rel_res| = %.3e\n', ...
    abs(max_body_residual - sol.rel_res));

if plots.make_figure
    residual_floor = get_positive_floor(body_residual);
    residual_log = log10(max(body_residual, residual_floor));
    default_residual_clim = expand_if_degenerate([min(residual_log), max(residual_log)]);

    fig = figure('Name','apr14 applied forcing','Color','w');
    ax = axes('Parent',fig);
    plot_disk_field(ax, q, geom.rad, body_forcing, plots.disk_points_full, ...
        plots.full_edge_color, plots.line_width_full);
    draw_domain_boundary(ax, geom_meta.L, geom.domain);
    cbar = colorbar(ax);
    apply_colormap_and_clim(ax, plots.force_colormap, plots.force_clim);
    style_colorbar(cbar, plots.font_size);
    ylabel(cbar, '$\sqrt{\|F\|_2^2 + (T/r)^2}$', ...
        'Interpreter', 'latex', 'FontSize', plots.font_size);
    xlabel(ax, '$x$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    ylabel(ax, '$y$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    style_axes(ax, plots.font_size);

    fig = figure('Name','apr14 maximum boundary speed','Color','w');
    ax = axes('Parent',fig);
    plot_disk_field(ax, q, geom.rad, body_speed, plots.disk_points_full, ...
        plots.full_edge_color, plots.line_width_full);
    draw_domain_boundary(ax, geom_meta.L, geom.domain);
    cbar = colorbar(ax);
    apply_colormap_and_clim(ax, plots.speed_colormap, plots.speed_clim);
    style_colorbar(cbar, plots.font_size);
    ylabel(cbar, '\textnormal{maximum boundary speed}', ...
        'Interpreter', 'latex', 'FontSize', plots.font_size);
    xlabel(ax, '$x$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    ylabel(ax, '$y$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    style_axes(ax, plots.font_size);

    fig = figure('Name','apr14 max boundary residual','Color','w');
    ax = axes('Parent',fig);
    plot_disk_field(ax, q, geom.rad, residual_log, plots.disk_points_full, ...
        plots.full_edge_color, plots.line_width_full);
    apply_colormap_and_clim(ax, plots.residual_colormap, ...
        choose_clim(plots.residual_clim, default_residual_clim));
    draw_domain_boundary(ax, geom_meta.L, geom.domain);
    cbar = colorbar(ax);
    style_colorbar(cbar, plots.font_size);
    ylabel(cbar, '$\log_{10}$ max relative boundary residual', ...
        'Interpreter', 'latex', 'FontSize', plots.font_size);  
    xlabel(ax, '$x$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    ylabel(ax, '$y$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    style_axes(ax, plots.font_size);
end

function value = final_gmres_residual(sol)
if isfield(sol,'resvec') && ~isempty(sol.resvec)
    value = sol.resvec(end);
else
    value = nan;
end
end

function [U, W] = unpack_UW(UW)
U = [UW(1:3:end), UW(2:3:end)];
W = UW(3:3:end);
end

function floor_value = get_positive_floor(values)
positive_values = values(values > 0);
if isempty(positive_values)
    floor_value = eps;
else
    floor_value = max(eps, 0.5*min(positive_values));
end
end

function limits = expand_if_degenerate(limits)
if limits(2) <= limits(1)
    limits = limits + [-0.5, 0.5];
end
end

function plot_disk_field(ax, q, rad, values, n_circle, edge_color, line_width)
theta = linspace(0, 2*pi, n_circle + 1);
theta(end) = [];
offsets = rad*(cos(theta) + 1i*sin(theta));
P = numel(q);

faces = reshape(1:P*n_circle, n_circle, P).';
vertices = zeros(P*n_circle, 2);
for k = 1:P
    idx = faces(k,:);
    zk = q(k) + offsets;
    vertices(idx,1) = real(zk);
    vertices(idx,2) = imag(zk);
end

patch(ax, 'Faces', faces, 'Vertices', vertices, ...
    'FaceVertexCData', values(:), ...
    'FaceColor', 'flat', ...
    'EdgeColor', edge_color, ...
    'LineWidth', line_width);

axis(ax, 'equal');
set(ax, 'Color', 'w');
box(ax, 'on');
end

function limits = choose_clim(user_limits, default_limits)
if isempty(user_limits)
    limits = default_limits;
else
    limits = user_limits;
end
end

function apply_colormap_and_clim(ax, cmap, clim_values)
colormap(ax, cmap);
if ~isempty(clim_values)
    caxis(ax, clim_values);
end
end

function style_colorbar(cbar, font_size)
cbar.TickLabelInterpreter = 'latex';
cbar.FontSize = font_size;
end

function style_axes(ax, font_size)
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', font_size);
grid(ax, 'on');
end

function draw_domain_boundary(ax, L, domain)
half_L = L/2;
style = 'k--';
if strcmp(domain,'periodic')
    style = 'r--';
end

hold(ax, 'on');
plot(ax, [-half_L, half_L, half_L, -half_L, -half_L], ...
    [-half_L, -half_L, half_L, half_L, -half_L], ...
    style, 'LineWidth', 1.0);
end
