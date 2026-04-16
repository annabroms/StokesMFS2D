clear;
close all;
clc;

% Compare mobility solves on random P-circle geometries using peanut
% compression and 1-body preconditioning.
% The residual figure shows the final GMRES estimated relative residual.
% The postprocessed surface residual is also stored in results.
% A third figure compares the number of GMRES unknowns.
% A fourth figure compares the ratio of 1-body unknowns to peanut unknowns.

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

script_name = mfilename;
script_date = 'Apr 14, 2026';

fprintf('=== Compare Cluster Iterations (Apr 14, 2026) ===\n');

rng_seed = 8;
rng(rng_seed);

P = 100;
rad = 1;
nruns = 5;
deltavec = logspace(-3,-1,10);

geom.generator = 'random_mc'; % 'cluster_gen' or 'random_mc'
geom.domain = 'boxed';          % used for random_mc
geom.phi = 0.65;                % used for random_mc
geom.n_sweeps = 30;             % used for random_mc
geom.show_generation_plot = false;
geom_tag = get_geometry_tag(geom);

N_c = 60; %80
N_f = 150;
N_peanut = 400;

io.run_experiment = 1;

opt_template = get2Dparams(P,N_c,N_f);
opt_template.rad = rad;
opt_template.delta_pair = 0.2;
opt_template.N_peanut = N_peanut;
opt_template.gmres_tol = 1e-8;
opt_template.maxit = 2000;
opt_template.visualise_sol = 0;
opt_template.visualise_grid = 0;
opt_template.debug = 0;
opt_template.gmres_verbose = 0;
opt_template.surface_error_mode = 'rel';
opt_template.reuse_pair_basis_by_sep = false;
opt_template.show_counter = 1;
opt_template.cmap = 1;
opt_template.self_correct = 1;
opt_template.use_dense = 1;
opt_template.get_bndry_field = 1;
opt_template.parallel_precomp = 1;

methods = struct( ...
    'name',{'peanut','one_body'}, ...
    'label',{'Peanut compression','1-body preconditioning'}, ...
    'solver',{@solve_mob_peanut_enhanced,@solve_mob_1B_enhanced}, ...
    'color',{[0.8500 0.3250 0.0980],[0 0.4470 0.7410]});

% Data I/O and plotting switches.
% Example 1: run and save without plotting
%   io.run_experiment = true;  io.save_results = true;
%   io.load_results = false;   plots.make_figures = false;
% Example 2: load saved results and plot without rerunning
%   io.run_experiment = false; io.save_results = false;
%   io.load_results = true;    plots.make_figures = true;
if io.run_experiment
    io.load_results = false;
    io.save_results = true;
    plots.make_figures = false;
else
    io.load_results = true;
    io.save_results = false;
    plots.make_figures = true;
end
io.data_filename = sprintf('%s_P%d_%s.mat', script_name, P, geom_tag);

if io.run_experiment && io.load_results
    error('apr14_compare_cluster_iterations:ConflictingIO', ...
        'Choose either io.run_experiment or io.load_results, not both.');
end
if ~io.run_experiment && ~io.load_results
    error('apr14_compare_cluster_iterations:NoInputSource', ...
        'Set either io.run_experiment or io.load_results to true.');
end

data_path = fullfile(repo_root,'data',io.data_filename);

if io.load_results
    if ~isfile(data_path)
        error('apr14_compare_cluster_iterations:MissingDataFile', ...
            'Could not find saved results file: %s', data_path);
    end
    loaded = load(data_path,'results');
    if ~isfield(loaded,'results')
        error('apr14_compare_cluster_iterations:BadDataFile', ...
            'The data file %s does not contain a ''results'' struct.', data_path);
    end
    results = loaded.results;
    fprintf('Loaded results from %s\n', data_path);
else
    ndelta = numel(deltavec);
    nmethods = numel(methods);

    iters = nan(ndelta,nruns,nmethods);
    gmres_residuals = nan(ndelta,nruns,nmethods);
    surface_residuals = nan(ndelta,nruns,nmethods);
    gmres_unknowns = nan(ndelta,nruns,nmethods);
    unknown_ratio = nan(ndelta,nruns);
    clusters = cell(ndelta,nruns);
    geom_meta = cell(ndelta,nruns);
    velocities = cell(ndelta,nruns);
    rotations = cell(ndelta,nruns);

    fprintf(['Running mobility comparison with geometry=%s, P=%d, nruns=%d, ', ...
        'N_c=%d, N_f=%d, N_peanut=%d\n'], ...
        geom.generator,P,nruns,N_c,N_f,N_peanut);

    for idelta = 1:ndelta
        delta = deltavec(idelta);
        fprintf('\nDelta %2d/%2d: %.3e\n', idelta, ndelta, delta);

        for irun = 1:nruns
            fprintf('  Run %2d/%2d\n', irun, nruns);

            [q, geom_meta_run] = generate_geometry_sample(P,delta,rad,geom);
            U = randn(P,2);
            W = randn(P,1);

            clusters{idelta,irun} = q;
            geom_meta{idelta,irun} = geom_meta_run;
            velocities{idelta,irun} = U;
            rotations{idelta,irun} = W;

            for imethod = 1:nmethods
                opt_run = opt_template;
                [~,sol] = methods(imethod).solver(q,U,W,opt_run);

                iters(idelta,irun,imethod) = sol.it;
                gmres_residuals(idelta,irun,imethod) = sol.resvec(end);
                surface_residuals(idelta,irun,imethod) = sol.rel_res;
                gmres_unknowns(idelta,irun,imethod) = sol.gmres_unknowns;

                fprintf('    %-22s it = %4d   unknowns = %6.0f   gmres = %.3e   surface = %.3e\n', ...
                    methods(imethod).label, sol.it, ...
                    gmres_unknowns(idelta,irun,imethod), ...
                    gmres_residuals(idelta,irun,imethod), ...
                    surface_residuals(idelta,irun,imethod));
            end

            unknown_ratio(idelta,irun) = gmres_unknowns(idelta,irun,2) / gmres_unknowns(idelta,irun,1);
        end
    end

    results = struct();
    results.script_name = script_name;
    results.script_date = script_date;
    results.rng_seed = rng_seed;
    results.P = P;
    results.rad = rad;
    results.nruns = nruns;
    results.geom = geom;
    results.geom_tag = geom_tag;
    results.deltavec = deltavec;
    results.method_names = {methods.name};
    results.method_labels = {methods.label};
    results.options = opt_template;
    results.iters = iters;
    results.gmres_residuals = gmres_residuals;
    results.surface_residuals = surface_residuals;
    results.gmres_unknowns = gmres_unknowns;
    results.unknown_ratio = unknown_ratio;
    results.clusters = clusters;
    results.geom_meta = geom_meta;
    results.velocities = velocities;
    results.rotations = rotations;

    if io.save_results
        save(data_path,'results');
        fprintf('\nSaved raw results to %s\n', data_path);
    end
end

results = ensureDerivedResults(results,methods);

if plots.make_figures
    plotClusterIterationResults(results,methods);
end

function [min_vals,max_vals,mean_vals] = summariseRuns(values)
nrows = size(values,1);
min_vals = nan(nrows,1);
max_vals = nan(nrows,1);
mean_vals = nan(nrows,1);

for i = 1:nrows
    row = values(i,:);
    row = row(isfinite(row));
    if isempty(row)
        continue
    end
    min_vals(i) = min(row);
    max_vals(i) = max(row);
    mean_vals(i) = mean(row);
end
end

function plotBand(x,lower,upper,color)
x = x(:).';
lower = lower(:).';
upper = upper(:).';

mask = isfinite(x) & isfinite(lower) & isfinite(upper);
if ~any(mask)
    return
end

x = x(mask);
lower = lower(mask);
upper = upper(mask);

fill([x, fliplr(x)], [upper, fliplr(lower)], color, ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');
end

function y = clampPositive(y)
mask = isfinite(y) & y > 0;
if ~any(mask)
    return
end
floor_val = min(y(mask));
y(~mask) = floor_val;
end

function tag = get_geometry_tag(geom)
generator = lower(char(string(geom.generator)));
switch generator
    case 'cluster_gen'
        tag = 'cluster_gen';
    case 'random_mc'
        tag = sprintf('random_mc_phi%.3f', geom.phi);
    otherwise
        error('apr14_compare_cluster_iterations:BadGenerator', ...
            'Unknown geom.generator = %s. Use ''cluster_gen'' or ''random_mc''.', ...
            geom.generator);
end
end

function [q, meta] = generate_geometry_sample(P,delta,rad,geom)
generator = lower(char(string(geom.generator)));

switch generator
    case 'cluster_gen'
        q = grow_cluster(P,delta,2,rad,[],geom.show_generation_plot,false);
        meta = build_cluster_gen_meta(q,rad,delta);

    case 'random_mc'
        geom_opt = struct();
        geom_opt.domain = geom.domain;
        geom_opt.phi = geom.phi;
        geom_opt.rad = rad;
        geom_opt.min_gap = delta;
        geom_opt.n_sweeps = geom.n_sweeps;
        geom_opt.visualise = geom.show_generation_plot;
        [q, meta] = random_discs_mc(P,geom_opt);
        meta.generator = 'random_mc';

    otherwise
        error('apr14_compare_cluster_iterations:BadGenerator', ...
            'Unknown geom.generator = %s. Use ''cluster_gen'' or ''random_mc''.', ...
            geom.generator);
end
end

function meta = build_cluster_gen_meta(q,rad,delta)
min_center_distance = get_min_center_distance_nonperiodic(q);

meta = struct();
meta.L = nan;
meta.phi = nan;
meta.phi_target = nan;
meta.phi_error = nan;
meta.phi_rel_error = nan;
meta.rad = rad;
meta.domain = 'cluster_gen';
meta.min_gap = delta;
meta.d_min_allowed = 2*rad + delta;
meta.n_sweeps = nan;
meta.acceptance_by_sweep = [];
meta.move_scale_by_sweep = [];
meta.min_center_distance = min_center_distance;
meta.min_surface_gap = min_center_distance - 2*rad;
meta.rng_seed = [];
meta.target_acceptance = nan;
meta.move_scale0 = nan;
meta.final_move_scale = nan;
meta.final_acceptance = nan;
meta.initial_layout = 'grow_cluster';
meta.initial_layout_details = struct();
meta.d_init_min = min_center_distance;
meta.generator = 'cluster_gen';
end

function dmin = get_min_center_distance_nonperiodic(q)
P = numel(q);
dmin = inf;
for j = 2:P
    dmin = min(dmin, min(abs(q(j) - q(1:j-1))));
end
end

function results = ensureDerivedResults(results,methods)
nmethods = numel(methods);

if ~isfield(results,'unknown_ratio') || isempty(results.unknown_ratio)
    results.unknown_ratio = results.gmres_unknowns(:,:,2) ./ results.gmres_unknowns(:,:,1);
end

summary = repmat(struct( ...
    'name','', ...
    'label','', ...
    'iter_min',[], ...
    'iter_max',[], ...
    'iter_mean',[], ...
    'gmres_min',[], ...
    'gmres_max',[], ...
    'gmres_mean',[], ...
    'unknowns_min',[], ...
    'unknowns_max',[], ...
    'unknowns_mean',[], ...
    'unknown_ratio_min',[], ...
    'unknown_ratio_max',[], ...
    'unknown_ratio_mean',[], ...
    'surface_min',[], ...
    'surface_max',[], ...
    'surface_mean',[]), nmethods, 1);

for imethod = 1:nmethods
    summary(imethod).name = methods(imethod).name;
    summary(imethod).label = methods(imethod).label;
    [summary(imethod).iter_min,summary(imethod).iter_max,summary(imethod).iter_mean] = ...
        summariseRuns(results.iters(:,:,imethod));
    [summary(imethod).gmres_min,summary(imethod).gmres_max,summary(imethod).gmres_mean] = ...
        summariseRuns(results.gmres_residuals(:,:,imethod));
    [summary(imethod).unknowns_min,summary(imethod).unknowns_max,summary(imethod).unknowns_mean] = ...
        summariseRuns(results.gmres_unknowns(:,:,imethod));
    [summary(imethod).surface_min,summary(imethod).surface_max,summary(imethod).surface_mean] = ...
        summariseRuns(results.surface_residuals(:,:,imethod));
end

[unknown_ratio_min,unknown_ratio_max,unknown_ratio_mean] = summariseRuns(results.unknown_ratio);
for imethod = 1:nmethods
    summary(imethod).unknown_ratio_min = unknown_ratio_min;
    summary(imethod).unknown_ratio_max = unknown_ratio_max;
    summary(imethod).unknown_ratio_mean = unknown_ratio_mean;
end

results.summary = summary;
end

function plotClusterIterationResults(results,methods)
deltavec = results.deltavec;
summary = results.summary;
nmethods = numel(methods);

figure(1);
clf;
hold on;
for imethod = 1:nmethods
    plotBand(deltavec,summary(imethod).iter_min,summary(imethod).iter_max,methods(imethod).color);
    plot(deltavec,summary(imethod).iter_mean,'-o', ...
        'Color',methods(imethod).color, ...
        'LineWidth',1.5, ...
        'MarkerFaceColor',methods(imethod).color, ...
        'DisplayName',methods(imethod).label);
end
set(gca,'XScale','log');
set(gca,'TickLabelInterpreter','latex');
grid on;
axis tight;
xlabel('$\delta$','Interpreter','latex');
ylabel('GMRES iterations','Interpreter','latex');
title('Cluster solve iterations','Interpreter','latex');
legend('Location','best','Interpreter','latex');

figure(2);
clf;
hold on;
for imethod = 1:nmethods
    gmres_min = clampPositive(summary(imethod).gmres_min);
    gmres_max = clampPositive(summary(imethod).gmres_max);
    gmres_mean = clampPositive(summary(imethod).gmres_mean);
    plotBand(deltavec,gmres_min,gmres_max,methods(imethod).color);
    plot(deltavec,gmres_mean,'-o', ...
        'Color',methods(imethod).color, ...
        'LineWidth',1.5, ...
        'MarkerFaceColor',methods(imethod).color, ...
        'DisplayName',methods(imethod).label);
end
set(gca,'XScale','log','YScale','log');
set(gca,'TickLabelInterpreter','latex');
grid on;
axis tight;
xlabel('$\delta$','Interpreter','latex');
ylabel('Final GMRES estimated relative residual','Interpreter','latex');
title('Cluster solve residuals','Interpreter','latex');
legend('Location','best','Interpreter','latex');

figure(3);
clf;
hold on;
for imethod = 1:nmethods
    plotBand(deltavec,summary(imethod).unknowns_min,summary(imethod).unknowns_max,methods(imethod).color);
    plot(deltavec,summary(imethod).unknowns_mean,'-o', ...
        'Color',methods(imethod).color, ...
        'LineWidth',1.5, ...
        'MarkerFaceColor',methods(imethod).color, ...
        'DisplayName',methods(imethod).label);
end
set(gca,'XScale','log');
set(gca,'TickLabelInterpreter','latex');
grid on;
axis tight;
xlabel('$\delta$','Interpreter','latex');
ylabel('GMRES unknowns','Interpreter','latex');
title('Cluster GMRES unknown counts','Interpreter','latex');
legend('Location','best','Interpreter','latex');

figure(4);
clf;
hold on;
plotBand(deltavec,summary(1).unknown_ratio_min,summary(1).unknown_ratio_max,methods(2).color);
plot(deltavec,summary(1).unknown_ratio_mean,'-o', ...
    'Color',methods(2).color, ...
    'LineWidth',1.5, ...
    'MarkerFaceColor',methods(2).color, ...
    'DisplayName','1-body / peanut');
set(gca,'XScale','log');
set(gca,'TickLabelInterpreter','latex');
grid on;
axis tight;
xlabel('$\delta$','Interpreter','latex');
ylabel('Unknown ratio','Interpreter','latex');
title('1-body to peanut GMRES unknown ratio','Interpreter','latex');
legend('Location','best','Interpreter','latex');
end
