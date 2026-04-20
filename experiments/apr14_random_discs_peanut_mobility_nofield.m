close all;
clc;

script_name = mfilename;
script_date = 'Apr 14, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Random Discs Peanut Mobility (Apr 14, 2026) ===\n');

% Geometry
geom.P = 250;
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
solver.get_precomp_time = true;

% Comparison figure
compare.make_triptych_figure = true;
compare.reference_N_c = 120;
compare.generate_reference_if_missing = true;
compare.figure_width_cm = 14.0;

% Data I/O
% Use io.save_results = true together with plots.make_figure = false to
% run the solve and store raw data without plotting. Later, set
% io.load_results = true to load the saved data and plot without rerunning.
run_test = 1; 
if run_test
    io.load_results = false;
    io.save_results = true;
    plots.make_figure = false;
else
    io.load_results = true;
    io.save_results = false;
    plots.make_figure = true;
end

io.data_filename = build_results_filename(script_name, geom.P, geom.phi, solver.N_c);
io.name2 = io.data_filename; %build_results_filename(script_name, geom.P, geom.phi, compare.reference_N_c);

% Reporting and visualisation
report.count_close_pairs = (geom.P <= 3000) && strcmp(geom.domain,'boxed');
plots.font_size = 16;
plots.full_edge_color = [0.10 0.10 0.10];
plots.line_width_full = 0.15;
plots.force_colormap = turbo(256);
cmap = winter(256);
cmap = cmap.^1.2;                    % darken
cmap = cmap + 0.15*(1 - cmap);       % brighten top end
cmap = min(cmap,1);
plots.speed_colormap = cmap;
%plots.speed_colormap = winter(256);
%plots.residual_colormap = hot(256);
cmap = copper(256);
% Stretch contrast
% alpha = 4;              % >1 = more contrast
% cmap = (cmap - 0.5)*alpha + 0.5;
% % Clip to valid range
% cmap = max(0, min(1, cmap));
cmap = cmap.^0.7;
plots.residual_colormap = cmap;
plots.force_clim = [];
plots.speed_clim = [];
plots.residual_clim = [];
plots.error_colormap = turbo(256);
plots.error_colormap = summer(256);
plots.error_clim = [];
plots.pair_gap_bins = 40;
plots.pair_gap_face_color = [0.20 0.20 0.24];
plots.pair_gap_edge_color = [0.05 0.05 0.05];

% TikZ export
export.write_triptych_tikz = true;
export.output_dir = fullfile(repo_root,'visuals','tikz');
export.max_particles_for_direct_tikz = 200;
export.colorbar_samples = 96;
export.large_P_resolution = 400;
export.basename = sprintf('%s_triptych_P%d_phi%.3f_Nc%d_vs_Nc%d', ...
    script_name, geom.P, geom.phi, solver.N_c, compare.reference_N_c);

if geom.P <= 300
    plots.disk_points_full = 56;
elseif geom.P <= 1500
    plots.disk_points_full = 36;
elseif geom.P <= 5000
    plots.disk_points_full = 24;
else
    plots.disk_points_full = 18;
end

if io.load_results && io.save_results
    error('apr14_random_discs_peanut_mobility:ConflictingIO', ...
        'Choose either io.load_results or io.save_results, not both.');
end

main_data_path = fullfile(repo_root,'data',io.data_filename);
reference_data_path = fullfile(repo_root,'data',io.name2);

if io.load_results
    if ~isfile(main_data_path)
        error('apr14_random_discs_peanut_mobility:MissingDataFile', ...
            'Could not find saved results file: %s', main_data_path);
    end
    results = load_results_file(main_data_path);
    fprintf('Loaded results from %s\n', main_data_path);
else
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

    opt = build_solver_opt(geom.P, solver, geom);
    opt.get_bndry_field = 0;
    opt.RAM_check = 0; 
    

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

    results = struct();
    results.script_name = script_name;
    results.script_date = script_date;
    results.geom = geom;
    results.mc = mc;
    results.loads = loads;
    results.solver = solver;
    results.report = report;
    results.geom_meta = geom_meta;
    results.q = q;
    results.n_close = n_close;
    results.F = F;
    results.T = T;
    results.UW = UW;
    results.sol = sol;
    results.solve_time = solve_time;
    results.body_forcing = body_forcing;
    results.body_speed = body_speed;
    results.body_residual = body_residual;
    results.gmres_residual = gmres_residual;
    results.opt = opt;

    if io.save_results
        save(main_data_path,'results');
        fprintf('Saved results to %s\n', main_data_path);
    end
end

need_reference_results = compare.make_triptych_figure || export.write_triptych_tikz;
reference_results = [];
body_velocity_error = [];
if need_reference_results
    if ~isfile(reference_data_path)
        if compare.generate_reference_if_missing
            reference_results = generate_reference_results(results, compare.reference_N_c);
            save_results_file(reference_data_path, reference_results);
        else
            error('apr14_random_discs_peanut_mobility:MissingReferenceFile', ...
                ['Could not find reference results file: %s\n' ...
                 'Generate it with N_c = %d or set compare.generate_reference_if_missing = true.'], ...
                reference_data_path, compare.reference_N_c);
        end
    end
    if isempty(reference_results)
        reference_results = load_results_file(reference_data_path);
    end
    validate_reference_results(results, reference_results, compare.reference_N_c);
    body_velocity_error = max_boundary_velocity_error(results.UW, reference_results.UW, results.geom.rad);
end

geom_run = results.geom;
solver_run = results.solver;
geom_meta = results.geom_meta;
q = results.q;
n_close = results.n_close;
F = results.F;
T = results.T;
UW = results.UW;
sol = results.sol;
solve_time = results.solve_time;
body_forcing = results.body_forcing;
body_speed = results.body_speed;
body_residual = results.body_residual;
gmres_residual = results.gmres_residual;
opt = results.opt;

fprintf('Geometry:\n');
fprintf('  domain=%s, P=%d, phi_target=%.6f, phi=%.6f, L=%.6f\n', ...
    geom_meta.domain, geom_run.P, geom_meta.phi_target, geom_meta.phi, geom_meta.L);
fprintf('  min allowed gap=%.3e, measured min gap=%.3e\n', ...
    geom_meta.min_gap, geom_meta.min_surface_gap);
if isfinite(n_close)
    fprintf('  close pairs below delta_pair=%.3f: %d\n', solver_run.delta_pair, n_close);
end

fprintf('Solver:\n');
fprintf('  N_c=%d, N_f=%d, N_peanut=%d, gmres_tol=%.1e\n', ...
    opt.N_c, opt.N_f, opt.N_peanut, opt.gmres_tol);
fprintf('  iterations=%d, unknowns=%d, gmres_res=%.3e, surf_rel=%.3e, time=%.2fs\n', ...
    sol.it, sol.gmres_unknowns, gmres_residual, sol.rel_res, solve_time);
if isfield(sol,'precomp_time') && isstruct(sol.precomp_time) && isfinite(sol.precomp_time.total)
    fprintf(['  precomp total=%.2fs, one-body=%.2fs, ', ...
        'two-body/peanut=%.2fs\n'], ...
        sol.precomp_time.total, ...
        sol.precomp_time.one_body, ...
        sol.precomp_time.two_body_or_peanut);
end

if ~plots.make_figure
    return;
end

if isempty(body_residual)
    error('apr14_random_discs_peanut_mobility:MissingResiduals', ...
        'solve_mob_peanut_enhanced did not return sol.body_rel_res_max.');
end

max_body_residual = max(body_residual);
fprintf('  max(body_rel_res_max)=%.3e\n', max_body_residual);
fprintf('  |max(body_rel_res_max) - sol.rel_res| = %.3e\n', ...
    abs(max_body_residual - sol.rel_res));
if ~isempty(body_velocity_error)
    fprintf('Comparison:\n');
    fprintf('  reference file: %s\n', io.name2);
    fprintf('  reference N_c=%d, max boundary velocity error=%.3e\n', ...
        compare.reference_N_c, max(body_velocity_error));
end

residual_floor = get_positive_floor(body_residual);
residual_log = log10(max(body_residual, residual_floor));
default_residual_clim = expand_if_degenerate([min(residual_log), max(residual_log)]);

triptych_panels = [];
triptych_fig = [];
if ~isempty(body_velocity_error)
    default_speed_clim = expand_if_degenerate([min(body_speed), max(body_speed)]);
    default_speed_clim(1) = 0;
    default_error_clim = expand_if_degenerate([min(body_velocity_error), max(body_velocity_error)]);
    error_floor = get_positive_floor(body_velocity_error);
    error_log = log10(max(body_velocity_error, error_floor));
    default_error_clim = expand_if_degenerate([min(error_log), max(error_log)]);
    triptych_panels = make_triptych_panels( ...
        body_speed, choose_clim(plots.speed_clim, default_speed_clim), plots.speed_colormap, ...
        '$\mathrm{max\ boundary\ speed}$', ...
        residual_log, choose_clim(plots.residual_clim, default_residual_clim), copper(256), ...
        '$\mathrm{max\ boundary\ rel\ res}$', ...
        error_log, choose_clim(plots.error_clim, default_error_clim), plasma(256), ...
        '$\mathrm{max\ boundary\ vel\ err}$');
end

if plots.make_figure

    fig = figure('Name','apr14 applied forcing','Color','w');
    ax = axes('Parent',fig);
    plot_disk_field(ax, q, geom_run.rad, body_forcing, plots.disk_points_full, ...
        plots.full_edge_color, plots.line_width_full);
    %draw_domain_boundary(ax, geom_meta.L, geom_run.domain);
    cbar = colorbar(ax);
    apply_colormap_and_clim(ax, plots.force_colormap, plots.force_clim);
    style_colorbar(cbar, plots.font_size);
    ylabel(cbar, '$\sqrt{\|F\|_2^2 + (T/r)^2}$', ...
        'Interpreter', 'latex', 'FontSize', plots.font_size);
    xlabel(ax, '$x$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    ylabel(ax, '$y$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    style_axes(ax, plots.font_size);
    axis tight
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])

    fig = figure('Name','apr14 maximum boundary speed','Color','w');
    ax = axes('Parent',fig);
    plot_disk_field(ax, q, geom_run.rad, body_speed, plots.disk_points_full, ...
        plots.full_edge_color, plots.line_width_full);
    %draw_domain_boundary(ax, geom_meta.L, geom_run.domain);
    cbar = colorbar(ax);
    apply_colormap_and_clim(ax, plots.speed_colormap, plots.speed_clim);
    style_colorbar(cbar, plots.font_size);
    ylabel(cbar, '\textnormal{maximum boundary speed}', ...
        'Interpreter', 'latex', 'FontSize', plots.font_size);
    xlabel(ax, '$x$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    ylabel(ax, '$y$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    style_axes(ax, plots.font_size);
    axis tight
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])

    fig = figure('Name','apr14 max boundary residual','Color','w');
    ax = axes('Parent',fig);
    plot_disk_field(ax, q, geom_run.rad, residual_log, plots.disk_points_full, ...
        plots.full_edge_color, plots.line_width_full);
    apply_colormap_and_clim(ax, plots.residual_colormap, ...
        choose_clim(plots.residual_clim, default_residual_clim));
    %draw_domain_boundary(ax, geom_meta.L, geom_run.domain);
    cbar = colorbar(ax);
    style_colorbar(cbar, plots.font_size);
    ylabel(cbar, '$\log_{10}$ max relative boundary residual', ...
        'Interpreter', 'latex', 'FontSize', plots.font_size);  
    xlabel(ax, '$x$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    ylabel(ax, '$y$', 'Interpreter', 'latex', 'FontSize', plots.font_size);
    style_axes(ax, plots.font_size);
    axis tight
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])

    fig = figure('Name','apr14 pair-separation histogram','Color','w');
    ax = axes('Parent',fig);
    plot_pair_gap_histogram(ax, q, geom_run.rad, geom_run.domain, geom_meta.L, ...
        solver_run.delta_pair, plots.pair_gap_bins, ...
        plots.pair_gap_face_color, plots.pair_gap_edge_color);
    xlabel(ax, '$\delta$', ...
        'Interpreter', 'latex', 'FontSize', plots.font_size);
    ylabel(ax, '\textnormal{pair count}', ...
        'Interpreter', 'latex', 'FontSize', plots.font_size);
    style_axes(ax, plots.font_size);
    set(ax, 'XScale', 'log', 'XMinorGrid', 'on');

    if compare.make_triptych_figure && ~isempty(triptych_panels)
        triptych_fig = create_triptych_figure(q, geom_run.rad, triptych_panels, plots, compare, 'on', true);
    end
end

if export.write_triptych_tikz && ~isempty(triptych_panels)
    if geom_run.P <= export.max_particles_for_direct_tikz
        write_triptych_tikz(q, geom_run.rad, triptych_panels, compare, export);
    else
        export_fig = create_triptych_figure(q, geom_run.rad, triptych_panels, plots, compare, 'off', false);
        write_triptych_wrapper_tikz(export_fig, q, geom_run.rad, triptych_panels, compare, export);
        close(export_fig);
    end
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
    'EdgeColor', 'none', ...
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

function plot_pair_gap_histogram(ax, q, rad, domain, L, delta_pair, n_bins, face_color, edge_color)
P = numel(q);
close_gaps = [];
for j = 2:P
    diffs = q(j) - q(1:j-1);
    if strcmp(domain, 'periodic')
        dx = real(diffs);
        dy = imag(diffs);
        dx = dx - L*round(dx/L);
        dy = dy - L*round(dy/L);
        distances = hypot(dx, dy);
    else
        distances = abs(diffs);
    end
    gaps = distances - 2*rad;
    gaps = gaps(gaps < delta_pair);
    if ~isempty(gaps)
        close_gaps = [close_gaps; gaps(:)]; %#ok<AGROW>
    end
end

set(ax, 'XScale', 'log');
box(ax, 'on');
set(ax, 'Color', 'w');

if isempty(close_gaps)
    lower_edge = max(eps, 0.1*delta_pair);
    upper_edge = max(delta_pair, 1.05*lower_edge);
    histogram(ax, 'BinEdges', logspace(log10(lower_edge), log10(upper_edge), n_bins + 1), ...
        'BinCounts', zeros(1, n_bins), ...
        'FaceColor', face_color, 'EdgeColor', edge_color, 'LineWidth', 0.6);
    xlim(ax, [lower_edge, upper_edge]);
    text(ax, sqrt(lower_edge*upper_edge), 0, ...
        '\textnormal{no pairs with surface gap $< \delta_{\mathrm{pair}}$}', ...
        'Interpreter', 'latex', 'FontSize', ax.FontSize, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    return
end

lower_edge = max(eps, min(close_gaps));
upper_edge = max(delta_pair, 1.05*lower_edge);
if upper_edge <= lower_edge
    upper_edge = 1.05*lower_edge;
end
edges = logspace(log10(lower_edge), log10(upper_edge), n_bins + 1);

histogram(ax, close_gaps, 'BinEdges', edges, ...
    'FaceColor', face_color, 'EdgeColor', edge_color, 'LineWidth', 0.6);
xlim(ax, [lower_edge, upper_edge]);
end

function filename = build_results_filename(script_name, P, phi, N_c)
filename = sprintf('%s_P%d_phi%.3f_Nc%d_results.mat', ...
    script_name, P, phi, N_c);
end

function results = load_results_file(data_path)
loaded = load(data_path,'results');
if ~isfield(loaded,'results')
    error('apr14_random_discs_peanut_mobility:BadDataFile', ...
        'The data file %s does not contain a ''results'' struct.', data_path);
end
results = loaded.results;
end

function save_results_file(data_path, results)
save(data_path,'results');
fprintf('Saved results to %s\n', data_path);
end

function opt = build_solver_opt(P, solver, geom)
opt = get2Dparams(P, solver.N_c, solver.N_f);
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
opt.get_precomp_time = solver.get_precomp_time;
opt.parallel_precomp = 1;
end

function reference_results = generate_reference_results(results, reference_N_c)
fprintf('Generating reference results with N_c=%d\n', reference_N_c);

solver_ref = results.solver;
solver_ref.N_c = reference_N_c;
opt_ref = build_solver_opt(results.geom.P, solver_ref, results.geom);
if isfield(results,'opt') && isfield(results.opt,'parallel_precomp')
    opt_ref.parallel_precomp = results.opt.parallel_precomp;
end

tic;
[UW_ref, sol_ref] = solve_mob_peanut_enhanced(results.q, results.F, results.T, opt_ref);
solve_time_ref = toc;

[U_ref, W_ref] = unpack_UW(UW_ref);
reference_results = struct();
reference_results.script_name = results.script_name;
reference_results.script_date = results.script_date;
reference_results.geom = results.geom;
reference_results.mc = results.mc;
reference_results.loads = results.loads;
reference_results.solver = solver_ref;
reference_results.report = results.report;
reference_results.geom_meta = results.geom_meta;
reference_results.q = results.q;
reference_results.n_close = results.n_close;
reference_results.F = results.F;
reference_results.T = results.T;
reference_results.UW = UW_ref;
reference_results.sol = sol_ref;
reference_results.solve_time = solve_time_ref;
reference_results.body_forcing = results.body_forcing;
reference_results.body_speed = sqrt(sum(U_ref.^2,2)) + results.geom.rad*abs(W_ref);
reference_results.body_residual = sol_ref.body_rel_res_max;
reference_results.gmres_residual = final_gmres_residual(sol_ref);
reference_results.opt = opt_ref;
end

function validate_reference_results(main_results, reference_results, reference_N_c)
tol = 1e-12;

if ~isfield(reference_results,'opt') || ~isfield(reference_results.opt,'N_c')
    error('apr14_random_discs_peanut_mobility:BadReferenceFile', ...
        'The reference file is missing opt.N_c.');
end

% if reference_results.opt.N_c ~= reference_N_c
%     error('apr14_random_discs_peanut_mobility:ReferenceNcMismatch', ...
%         'Expected reference N_c = %d, found %d.', ...
%         reference_N_c, reference_results.opt.N_c);
% end

if numel(main_results.q) ~= numel(reference_results.q) || ...
        max(abs(main_results.q - reference_results.q)) > tol
    error('apr14_random_discs_peanut_mobility:ReferenceGeometryMismatch', ...
        'The reference file does not use the same particle centers q.');
end

if max(abs(main_results.F(:) - reference_results.F(:))) > tol || ...
        max(abs(main_results.T(:) - reference_results.T(:))) > tol
    error('apr14_random_discs_peanut_mobility:ReferenceLoadMismatch', ...
        'The reference file does not use the same forcing/torque data.');
end
end

function body_error = max_boundary_velocity_error(UW, UW_ref, rad)
[U, W] = unpack_UW(UW);
[U_ref, W_ref] = unpack_UW(UW_ref);
body_error = sqrt(sum((U - U_ref).^2,2)) + rad*abs(W - W_ref);
end

function panels = make_triptych_panels( ...
    values1, clim1, cmap1, label1, ...
    values2, clim2, cmap2, label2, ...
    values3, clim3, cmap3, label3)
ticks1 = [clim1(1), clim1(2)];
labels1 = {'0', sprintf('%.1f', clim1(2))};
ticks2 = [clim2(1), clim2(2)];
labels2 = format_triptych_ticklabels(ticks2);
ticks3 = [clim3(1), clim3(2)];
labels3 = format_triptych_ticklabels(ticks3);

panels = struct( ...
    'values', {values1(:), values2(:), values3(:)}, ...
    'clim', {clim1, clim2, clim3}, ...
    'cmap', {cmap1, cmap2, cmap3}, ...
    'colorbar_label', {label1, label2, label3}, ...
    'tick_values', {ticks1, ticks2, ticks3}, ...
    'tick_labels', {labels1, labels2, labels3});
end

function fig = create_triptych_figure(q, rad, panels, plots, compare, visible_state, show_colorbar_text)
if nargin < 6 || isempty(visible_state)
    visible_state = 'on';
end
if nargin < 7
    show_colorbar_text = true;
end

global_bounds = compute_square_bounds(q, rad);
layout = get_triptych_layout(compare.figure_width_cm, global_bounds);
fig = figure('Name','apr14 speed-residual-error triptych', ...
    'Color','w', 'Visible', visible_state, 'Units','centimeters', ...
    'Position',[1, 1, compare.figure_width_cm, layout.figure_height]);

for k = 1:3
    slice_bounds = get_triptych_slice_bounds(global_bounds, k);
    ax = axes('Parent', fig, 'Units','centimeters', ...
        'Position', [layout.lefts(k), layout.axes_bottom, ...
                     layout.panel_width, layout.panel_height]);
    plot_disk_field(ax, q, rad, panels(k).values, plots.disk_points_full, ...
        plots.full_edge_color, plots.line_width_full);
    apply_colormap_and_clim(ax, panels(k).cmap, panels(k).clim);
    xlim(ax, slice_bounds(1:2));
    ylim(ax, global_bounds(3:4));
    style_triptych_axes(ax);

    cbar = colorbar(ax, 'southoutside');
    if show_colorbar_text
        style_colorbar(cbar, max(plots.font_size - 3, 8));
        cbar.Ticks = panels(k).tick_values;
        cbar.TickLabels = panels(k).tick_labels;
        cbar.TickDirection = 'in';
        xlabel(cbar, panels(k).colorbar_label, 'Interpreter', 'latex', ...
            'FontSize', max(plots.font_size - 3, 8));
    else
        cbar.Ticks = [];
        cbar.TickLabels = {};
        cbar.Label.String = '';
    end
    drawnow;
    cbar.Units = 'centimeters';
    cbar.Position = [layout.lefts(k) + layout.cb_side_gap, layout.cb_bottom, ...
        layout.panel_width - 2*layout.cb_side_gap, layout.cb_height];
    ax.Position = [layout.lefts(k), layout.axes_bottom, ...
        layout.panel_width, layout.panel_height];
end
end

function style_triptych_axes(ax)
axis(ax, 'equal');
box(ax, 'on');
grid(ax, 'off');
set(ax, 'XTick', [], 'YTick', [], ...
    'TickLabelInterpreter', 'latex', ...
    'XColor', 'k', 'YColor', 'k', ...
    'LineWidth', 0.6, ...
    'Color', 'w');
end

function layout = get_triptych_layout(fig_width_cm, global_bounds)
margin_x = 0.45;
gap_x = 0.0;
top_margin = 0.30;
cb_bottom = 0.95;
cb_height = 0.25;
cb_gap = 0.70;
cb_side_gap = 0.4; %0.18;
axes_bottom = cb_bottom + cb_height + cb_gap;
usable_width = fig_width_cm - 2*margin_x;
panel_width = usable_width/3;
slice_width = (global_bounds(2) - global_bounds(1))/3;
global_height = global_bounds(4) - global_bounds(3);
panel_height = panel_width * global_height / slice_width;
left0 = margin_x;
figure_height = axes_bottom + panel_height + top_margin;

layout = struct();
layout.figure_height = figure_height;
layout.panel_width = panel_width;
layout.panel_height = panel_height;
layout.axes_bottom = axes_bottom;
layout.cb_bottom = cb_bottom;
layout.cb_height = cb_height;
layout.cb_side_gap = cb_side_gap;
layout.lefts = left0 + (0:2)*(panel_width + gap_x);
end

function slice_bounds = get_triptych_slice_bounds(global_bounds, panel_index)
slice_width = (global_bounds(2) - global_bounds(1))/3;
xmin = global_bounds(1) + (panel_index - 1)*slice_width;
xmax = xmin + slice_width;
slice_bounds = [xmin, xmax, global_bounds(3), global_bounds(4)];
end

function write_triptych_tikz(q, rad, panels, compare, export)
if ~exist(export.output_dir,'dir')
    mkdir(export.output_dir);
end

global_bounds = compute_square_bounds(q, rad);
layout = get_triptych_layout(compare.figure_width_cm, global_bounds);
picture_file = fullfile(export.output_dir, [export.basename '_picture.tex']);
demo_file = fullfile(export.output_dir, [export.basename '_demo.tex']);

fid = fopen(picture_file,'w');
if fid < 0
    error('apr14_random_discs_peanut_mobility:TikzOpenFailed', ...
        'Could not open %s for writing.', picture_file);
end
cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'%% Auto-generated from apr14_random_discs_peanut_mobility.m\n');
fprintf(fid,'%% Adjust these font commands to match your paper.\n');
fprintf(fid,'\\providecommand{\\aprFourteenTriptychColorbarLabelFont}{\\fontsize{8}{9}\\selectfont}\n');
fprintf(fid,'\\providecommand{\\aprFourteenTriptychTickFont}{\\fontsize{7}{8}\\selectfont}\n');
fprintf(fid,'\\begin{tikzpicture}[x=1cm,y=1cm]\n');

for k = 1:numel(panels)
    slice_bounds = get_triptych_slice_bounds(global_bounds, k);
    write_triptych_panel_tikz(fid, q, rad, global_bounds, slice_bounds, panels(k), ...
        layout.lefts(k), layout.axes_bottom, layout.panel_width, layout.panel_height, ...
        layout.cb_bottom, layout.cb_height, export.colorbar_samples);
end

fprintf(fid,'\\end{tikzpicture}\n');

write_triptych_demo_tex(demo_file, [export.basename '_picture.tex']);
fprintf('Wrote TikZ picture: %s\n', picture_file);
fprintf('Wrote TikZ demo: %s\n', demo_file);
end

function write_triptych_wrapper_tikz(fig, q, rad, panels, compare, export)
if ~exist(export.output_dir,'dir')
    mkdir(export.output_dir);
end

global_bounds = compute_square_bounds(q, rad);
layout = get_triptych_layout(compare.figure_width_cm, global_bounds);
pdf_file = fullfile(export.output_dir, [export.basename '_image.pdf']);
picture_file = fullfile(export.output_dir, [export.basename '_picture.tex']);
demo_file = fullfile(export.output_dir, [export.basename '_demo.tex']);

exportgraphics(fig, pdf_file, 'ContentType', 'image', ...
    'Resolution', export.large_P_resolution);

fid = fopen(picture_file,'w');
if fid < 0
    error('apr14_random_discs_peanut_mobility:TikzOpenFailed', ...
        'Could not open %s for writing.', picture_file);
end
cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'%% Auto-generated from apr14_random_discs_peanut_mobility.m\n');
fprintf(fid,'%% Large-P export: the particle geometry is embedded as a PDF image,\n');
fprintf(fid,'%% while colorbar ticks and labels remain editable in TikZ.\n');
fprintf(fid,'\\providecommand{\\aprFourteenTriptychColorbarLabelFont}{\\fontsize{8}{9}\\selectfont}\n');
fprintf(fid,'\\providecommand{\\aprFourteenTriptychTickFont}{\\fontsize{7}{8}\\selectfont}\n');
fprintf(fid,'\\begin{tikzpicture}\n');
fprintf(fid,'  \\node[inner sep=0pt,outer sep=0pt] (img) at (0,0) {%%\n');
fprintf(fid,'    \\includegraphics[width=%.6fcm]{%s}%%\n', ...
    compare.figure_width_cm, [export.basename '_image.pdf']);
fprintf(fid,'  };\n');
fprintf(fid,'  \\begin{scope}[shift={(img.south west)},\n');
fprintf(fid,'    x={($(img.south east)-(img.south west)$)},\n');
fprintf(fid,'    y={($(img.north west)-(img.south west)$)}]\n');

for k = 1:numel(panels)
    write_triptych_overlay_tikz(fid, panels(k), layout, k);
end

fprintf(fid,'  \\end{scope}\n');
fprintf(fid,'\\end{tikzpicture}\n');

write_triptych_demo_tex(demo_file, [export.basename '_picture.tex']);
fprintf('Wrote TikZ wrapper picture: %s\n', picture_file);
fprintf('Wrote TikZ demo: %s\n', demo_file);
fprintf('Wrote triptych image PDF: %s\n', pdf_file);
end

function write_triptych_panel_tikz(fid, q, rad, global_bounds, slice_bounds, panel, left, bottom, panel_width, panel_height, cb_bottom, cb_height, n_colorbar)
rgb = map_values_to_rgb(panel.values, panel.cmap, panel.clim);
radius_tex = panel_width*rad/(slice_bounds(2) - slice_bounds(1));
fprintf(fid,'  \\fill[white] (%.6f, %.6f) rectangle (%.6f, %.6f);\n', ...
    left, bottom, left + panel_width, bottom + panel_height);
fprintf(fid,'  \\begin{scope}\n');
fprintf(fid,'    \\clip (%.6f, %.6f) rectangle (%.6f, %.6f);\n', ...
    left, bottom, left + panel_width, bottom + panel_height);

for j = 1:numel(q)
    x = left + panel_width*(real(q(j)) - slice_bounds(1))/(slice_bounds(2) - slice_bounds(1));
    y = bottom + panel_height*(imag(q(j)) - global_bounds(3))/(global_bounds(4) - global_bounds(3));
    fprintf(fid,['  \\fill[draw=none,fill=%s] (%.6f, %.6f) circle[radius=%.6fcm];\n'], ...
        tikz_rgb(rgb(j,:)), x, y, radius_tex);
end
fprintf(fid,'  \\end{scope}\n');
fprintf(fid,'  \\draw[line width=0.25pt] (%.6f, %.6f) rectangle (%.6f, %.6f);\n', ...
    left, bottom, left + panel_width, bottom + panel_height);

write_colorbar_tikz(fid, panel, left, cb_bottom, panel_width, cb_height, n_colorbar);
end

function write_triptych_overlay_tikz(fid, panel, layout, panel_index)
fig_width = layout.lefts(end) + layout.panel_width + layout.lefts(1);
fig_height = layout.figure_height;
cb_side_gap = layout.cb_side_gap;
cb_left = layout.lefts(panel_index) + cb_side_gap;
cb_width = layout.panel_width - 2*cb_side_gap;
cb_bottom = layout.cb_bottom;
cb_height = layout.cb_height;

for j = 1:numel(panel.tick_values)
    x = cb_left + cb_width*(panel.tick_values(j) - panel.clim(1))/(panel.clim(2) - panel.clim(1));
    fprintf(fid,'    \\draw[line width=0.2pt] (%.16f, %.16f) -- (%.16f, %.16f);\n', ...
        x/fig_width, cb_bottom/fig_height, x/fig_width, (cb_bottom - 0.08)/fig_height);
    fprintf(fid,'    \\node[font=\\aprFourteenTriptychTickFont,anchor=north] at (%.16f, %.16f) {%s};\n', ...
        x/fig_width, (cb_bottom - 0.10)/fig_height, panel.tick_labels{j});
end

fprintf(fid,'    \\node[font=\\aprFourteenTriptychColorbarLabelFont,anchor=north] at (%.16f, %.16f) {%s};\n', ...
    (cb_left + 0.5*cb_width)/fig_width, (cb_bottom - 0.42)/fig_height, panel.colorbar_label);
end

function write_colorbar_tikz(fid, panel, left, cb_bottom, cb_width, cb_height, n_colorbar)
cb_side_gap = 1; %0.18;
cb_left = left + cb_side_gap;
cb_width = cb_width - 2*cb_side_gap;
sample_values = linspace(panel.clim(1), panel.clim(2), n_colorbar + 1);
mid_values = 0.5*(sample_values(1:end-1) + sample_values(2:end));
rgb = map_values_to_rgb(mid_values(:), panel.cmap, panel.clim);

for j = 1:n_colorbar
    x0 = cb_left + cb_width*(j-1)/n_colorbar;
    x1 = cb_left + cb_width*j/n_colorbar;
    fprintf(fid,'  \\fill[draw=none,fill=%s] (%.6f, %.6f) rectangle (%.6f, %.6f);\n', ...
        tikz_rgb(rgb(j,:)), x0, cb_bottom, x1, cb_bottom + cb_height);
end

fprintf(fid,'  \\draw[line width=0.2pt] (%.6f, %.6f) rectangle (%.6f, %.6f);\n', ...
    cb_left, cb_bottom, cb_left + cb_width, cb_bottom + cb_height);

for j = 1:numel(panel.tick_values)
    x = cb_left + cb_width*(panel.tick_values(j) - panel.clim(1))/(panel.clim(2) - panel.clim(1));
    fprintf(fid,'  \\draw[line width=0.2pt] (%.6f, %.6f) -- (%.6f, %.6f);\n', ...
        x, cb_bottom, x, cb_bottom - 0.08);
    fprintf(fid,'  \\node[font=\\aprFourteenTriptychTickFont,anchor=north] at (%.6f, %.6f) {%s};\n', ...
        x, cb_bottom - 0.10, panel.tick_labels{j});
end

fprintf(fid,'  \\node[font=\\aprFourteenTriptychColorbarLabelFont,anchor=north] at (%.6f, %.6f) {%s};\n', ...
    cb_left + 0.5*cb_width, cb_bottom - 0.42, panel.colorbar_label);
end

function bounds = compute_square_bounds(q, rad)
xmin = min(real(q)) - rad;
xmax = max(real(q)) + rad;
ymin = min(imag(q)) - rad;
ymax = max(imag(q)) + rad;
center_x = 0.5*(xmin + xmax);
center_y = 0.5*(ymin + ymax);
half_width = 0.5*max(xmax - xmin, ymax - ymin);
bounds = [center_x - half_width, center_x + half_width, ...
          center_y - half_width, center_y + half_width];
end

function rgb = map_values_to_rgb(values, cmap, clim)
if clim(2) <= clim(1)
    clim = clim + [-0.5, 0.5];
end

t = (values(:) - clim(1)) / (clim(2) - clim(1));
t = min(max(t, 0), 1);
idx = 1 + t*(size(cmap,1) - 1);
i0 = floor(idx);
i1 = ceil(idx);
alpha = idx - i0;
rgb = (1 - alpha).*cmap(i0,:) + alpha.*cmap(i1,:);
end

function spec = tikz_rgb(rgb)
spec = sprintf('{rgb,1:red,%.6f;green,%.6f;blue,%.6f}', ...
    rgb(1), rgb(2), rgb(3));
end

function labels = format_triptych_ticklabels(values)
labels = cell(size(values));
for k = 1:numel(values)
    labels{k} = sprintf('%.1f', values(k));
end
end

function label = format_tex_number(value)
if value == 0
    label = '$0$';
    return
end

abs_value = abs(value);
if abs_value >= 1e-2 && abs_value < 1e2
    label = sprintf('$%.3g$', value);
else
    exponent = floor(log10(abs_value));
    mantissa = value / 10^exponent;
    label = sprintf('$%.2f\\times 10^{%d}$', mantissa, exponent);
end
end

function write_triptych_demo_tex(filename, picture_name)
fid = fopen(filename,'w');
if fid < 0
    error('apr14_random_discs_peanut_mobility:TikzDemoOpenFailed', ...
        'Could not open %s for writing.', filename);
end
cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'\\documentclass[tikz,border=4pt]{standalone}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
fprintf(fid,'\\usepackage{graphicx}\n');
fprintf(fid,'\\usepackage{tikz}\n');
fprintf(fid,'\\usetikzlibrary{calc}\n');
fprintf(fid,'\\begin{document}\n');
fprintf(fid,'\\input{%s}\n', picture_name);
fprintf(fid,'\\end{document}\n');
end

function cmap = plasma(n)
%PLASMA  Perceptually uniform purple–red–yellow colormap
%   cmap = plasma(n) returns an n-by-3 colormap.
%   If n is omitted, defaults to current figure colormap size.

if nargin < 1
    n = size(get(gcf,'colormap'),1);
end

% Key colors sampled from the original plasma colormap
base = [
    0.050383, 0.029803, 0.527975
    0.293478, 0.010213, 0.629490
    0.507860, 0.016333, 0.656202
    0.716387, 0.214982, 0.475290
    0.901807, 0.425087, 0.359688
    0.993033, 0.683153, 0.195911
    0.940015, 0.975158, 0.131326
];

% Interpolate to n colors
x = linspace(0,1,size(base,1));
xi = linspace(0,1,n);
cmap = interp1(x, base, xi, 'pchip');
end
