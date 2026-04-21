close all;
clc;

script_name = mfilename;
script_date = 'Apr 20, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')))
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== P Sweep Peanut Mobility (Apr 20, 2026) ===\n');

% Test parameters
phi = 0.65;
P_range = 50:50:1000;
%P_range = 5:5:20;
n_repeats = 10;
N_c = 60;
get_bndry_vel = true;

run_test = 1;

% Storage for results
n_p = length(P_range);
iterations_all = zeros(n_p, n_repeats);
residuals_all = zeros(n_p, n_repeats);
n_close_all = zeros(n_p, n_repeats);

% Create data directories
data_root = fullfile(repo_root, 'data');
individual_data_dir = fullfile(data_root, 'random_close_packing');
aggregated_filename = [script_name '_force_results.mat'];
%aggregated_filename = [script_name '_results.mat'];
aggregated_filepath = fullfile(data_root, aggregated_filename);

if run_test
    if ~exist(data_root, 'dir')
        mkdir(data_root);
    end
    if ~exist(individual_data_dir, 'dir')
        mkdir(individual_data_dir);
    end

    % Run experiments
    fprintf('Running experiments...\n');
    for i = 1:n_p
        P = P_range(i);
        fprintf('  P = %u (%d/%d)\n', P, i, n_p);

        for j = 1:n_repeats
            [iter, res, n_close, run_data] = run_single_experiment(P, phi, N_c, get_bndry_vel);
            iterations_all(i, j) = iter;
            residuals_all(i, j) = res;
            n_close_all(i, j) = n_close;

            % Save individual run data
            individual_filename = sprintf('P_%u_phi_%.3f_repeat_%d_results.mat', P, phi, j);
            individual_filepath = fullfile(individual_data_dir, individual_filename);
            save(individual_filepath, 'run_data');
        end
    end

    % Compute statistics
    iterations_mean = mean(iterations_all, 2);
    iterations_min = min(iterations_all, [], 2);
    iterations_max = max(iterations_all, [], 2);

    residuals_mean = mean(residuals_all, 2);
    residuals_min = min(residuals_all, [], 2);
    residuals_max = max(residuals_all, [], 2);

    n_close_mean = mean(n_close_all, 2);
    n_close_min = min(n_close_all, [], 2);
    n_close_max = max(n_close_all, [], 2);

    % Save aggregated results
    aggregated_results = struct();
    aggregated_results.script_name = script_name;
    aggregated_results.script_date = script_date;
    aggregated_results.P_range = P_range;
    aggregated_results.N_c = N_c;
    aggregated_results.get_bndry_vel = get_bndry_vel;
    aggregated_results.phi = phi;
    aggregated_results.n_repeats = n_repeats;
    aggregated_results.iterations_all = iterations_all;
    aggregated_results.residuals_all = residuals_all;
    aggregated_results.n_close_all = n_close_all;
    aggregated_results.iterations_mean = iterations_mean;
    aggregated_results.iterations_min = iterations_min;
    aggregated_results.iterations_max = iterations_max;
    aggregated_results.residuals_mean = residuals_mean;
    aggregated_results.residuals_min = residuals_min;
    aggregated_results.residuals_max = residuals_max;
    aggregated_results.n_close_mean = n_close_mean;
    aggregated_results.n_close_min = n_close_min;
    aggregated_results.n_close_max = n_close_max;

    save(aggregated_filepath, 'aggregated_results');
    fprintf('Saved aggregated results to %s\n', aggregated_filepath);

else
    if ~isfile(aggregated_filepath)
        fprintf('Aggregated results file not found. Extracting from individual run files...\n');
        iterations_all = extract_iterations_from_files(data_root, P_range, n_repeats);
        residuals_all = extract_residuals_from_files(data_root, P_range, n_repeats);
        n_close_all = extract_n_close_from_files(data_root, P_range, n_repeats);
    else
        f = load(aggregated_filepath);
        iterations_all = f.aggregated_results.iterations_all;

        % Try to load residuals_all, otherwise recompute from individual files
        if isfield(f.aggregated_results, 'residuals_all')
            residuals_all = f.aggregated_results.residuals_all;
        else
            residuals_all = extract_residuals_from_files(data_root, P_range, n_repeats);
        end

        % Try to load n_close_all, otherwise recompute from individual files
        if isfield(f.aggregated_results, 'n_close_all')
            n_close_all = f.aggregated_results.n_close_all;
        else
            % Fallback: recompute n_close from individual saved data
            n_close_all = extract_n_close_from_files(data_root, P_range, n_repeats);
        end
    end

    % Compute statistics
    iterations_mean = mean(iterations_all, 2);
    iterations_min = min(iterations_all, [], 2);
    iterations_max = max(iterations_all, [], 2);

    residuals_mean = mean(residuals_all, 2);
    residuals_min = min(residuals_all, [], 2);
    residuals_max = max(residuals_all, [], 2);

    n_close_mean = mean(n_close_all, 2);
    n_close_min = min(n_close_all, [], 2);
    n_close_max = max(n_close_all, [], 2);
end

%%
% Figure 1: Iterations vs p
fig1 = figure('Name', 'Iterations vs P', 'Color', 'w');
ax1 = axes('Parent', fig1);
hold(ax1, 'on');
plot(ax1, P_range, iterations_mean, 'b-o', 'LineWidth', 2, 'DisplayName', 'Mean');
%plot(ax1, P_range, iterations_min, 'g--s', 'LineWidth', 1.5, 'DisplayName', 'Min');
%plot(ax1, P_range, iterations_max, 'r--^', 'LineWidth', 1.5, 'DisplayName', 'Max');
fill_between_wrapper(ax1, P_range, iterations_min, iterations_max, 0.15);
hold(ax1, 'off');
xlabel(ax1, '$P$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel(ax1, 'Number of iterations', 'Interpreter', 'latex', 'FontSize', 14);
%legend(ax1, 'Location', 'best');
grid(ax1, 'on');
set(ax1, 'TickLabelInterpreter', 'latex');
axis tight

% Figure 2: Residuals vs p
fig2 = figure('Name', 'Residuals vs P', 'Color', 'w');
ax2 = axes('Parent', fig2);
hold(ax2, 'on');
semilogy(ax2, P_range, residuals_mean, 'b-o', 'LineWidth', 2, 'DisplayName', 'Mean');
%semilogy(ax2, phi_range, residuals_min, 'g--s', 'LineWidth', 1.5, 'DisplayName', 'Min');
%semilogy(ax2, phi_range, residuals_max, 'r--^', 'LineWidth', 1.5, 'DisplayName', 'Max');
fill_between_log_wrapper(ax2, P_range, residuals_min, residuals_max, 0.15);
hold(ax2, 'off');
xlabel(ax2, '$P$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel(ax2, 'Max relative residual', 'Interpreter', 'latex', 'FontSize', 14);
%legend(ax2, 'Location', 'best');
grid(ax2, 'on');
set(ax2, 'TickLabelInterpreter', 'latex');
axis tight

% Figure 3: Number of near pairs vs P
fig3 = figure('Name', 'Near pairs vs P', 'Color', 'w');
ax3 = axes('Parent', fig3);
hold(ax3, 'on');
plot(ax3, P_range, n_close_mean./P_range', 'b-o', 'LineWidth', 2, 'DisplayName', 'Mean');
fill_between_wrapper(ax3, P_range, n_close_min./P_range', n_close_max./P_range', 0.15);
hold(ax3, 'off');
xlabel(ax3, '$P$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel(ax3, 'Number of near pairs per particle', 'Interpreter', 'latex', 'FontSize', 14);
%legend(ax3, 'Location', 'best');
grid(ax3, 'on');
set(ax3, 'TickLabelInterpreter', 'latex');
axis tight

fprintf('Experiments complete.\n');

%% Helper function
function [iter, res, n_close, run_data] = run_single_experiment(P, phi, N_c, get_bndry_vel)
% RUN_SINGLE_EXPERIMENT Run a single peanut mobility experiment
%   [iter, res, n_close, run_data] = run_single_experiment(P, phi, N_c, get_bndry_vel)
%
%   Inputs:
%     P           - Number of particles
%     phi         - Target packing fraction
%     N_c         - Number of collocation points per particle
%     get_bndry_vel - Boolean flag to compute boundary velocity field
%
%   Outputs:
%     iter        - Number of GMRES iterations
%     res         - Maximum relative residual
%     n_close     - Number of close pairs
%     run_data    - Struct containing all results from this run

repo_root = fileparts(fileparts(mfilename('fullpath')));

% Generate geometry directly
geom_opt = struct();
geom_opt.domain = 'boxed';
geom_opt.phi = phi;
geom_opt.rad = 1;
geom_opt.min_gap = 1e-3;
geom_opt.n_sweeps = 200;
geom_opt.rng_seed = randi([1, 10000], 1, 1);
geom_opt.visualise = false;

[q, ~] = random_discs_mc(P, geom_opt);

% Count close pairs
delta_pair = 0.2;
n_close = count_close_pairs(q, delta_pair, 1);

% Build solver options
opt = get2Dparams(P, N_c, 150);
opt.rad = 1;
opt.delta_pair = delta_pair;
opt.N_peanut = 400;
opt.gmres_tol = 1e-8;
opt.maxit = 1000;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.debug = 0;
opt.gmres_verbose = 0;
opt.surface_error_mode = 'rel';
opt.reuse_pair_basis_by_sep = false;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.get_bndry_field = double(get_bndry_vel);
opt.get_precomp_time = 0;
opt.parallel_precomp = 0;
opt.RAM_check = 0;

% Generate random loads and solve
rng(randi([1, 10000], 1, 1), 'twister');
% F = zeros(P, 2);
% T = ones(P, 1);
F = rand(P,2);
T = rand(P,1);
[UW, sol] = solve_mob_peanut_enhanced(q, F, T, opt);

% Extract results
iter = sol.it;
if isfield(sol, 'body_rel_res_max') && ~isempty(sol.body_rel_res_max)
    res = max(sol.body_rel_res_max);
    res = sol.abs_res;
else
    res = sol.rel_res;

end

% Package run data
run_data = struct();
run_data.P = P;
run_data.phi = phi;
run_data.N_c = N_c;
run_data.get_bndry_vel = get_bndry_vel;
run_data.q = q;
run_data.F = F;
run_data.T = T;
run_data.UW = UW;
run_data.sol = sol;
run_data.iterations = iter;
run_data.residual = res;
run_data.n_close = n_close;
run_data.opt = opt;

end

%% Helper function to extract iterations from individual files
function iterations_all = extract_iterations_from_files(data_root, P_range, n_repeats)
% Extract GMRES iterations from individual saved run data files
% Used as fallback when aggregated results file is missing
n_p = length(P_range);
iterations_all = zeros(n_p, n_repeats);
individual_data_dir = fullfile(data_root, 'random_close_packing');
phi = 0.65;  % From script parameters

for i = 1:n_p
    P = P_range(i);
    for j = 1:n_repeats
        individual_filename = sprintf('P_%u_phi_%.3f_repeat_%d_results.mat', P, phi, j);
        individual_filepath = fullfile(individual_data_dir, individual_filename);

        if isfile(individual_filepath)
            loaded = load(individual_filepath, 'run_data');
            if isfield(loaded.run_data, 'iterations') && ~isempty(loaded.run_data.iterations)
                iterations_all(i, j) = loaded.run_data.iterations;
            elseif isfield(loaded.run_data, 'sol')
                % Fallback: extract from sol struct
                sol = loaded.run_data.sol;
                if isfield(sol, 'it')
                    iterations_all(i, j) = sol.it;
                end
            end
        else
            warning('Could not find file: %s', individual_filepath);
        end
    end
end
end

%% Helper function to extract residuals from individual files
function residuals_all = extract_residuals_from_files(data_root, P_range, n_repeats)
% Extract relative residuals from individual saved run data files
% Used as fallback when residuals were not stored in aggregated results
n_p = length(P_range);
residuals_all = zeros(n_p, n_repeats);
individual_data_dir = fullfile(data_root, 'random_close_packing');
phi = 0.65;  % From script parameters

for i = 1:n_p
    P = P_range(i);
    for j = 1:n_repeats
        individual_filename = sprintf('P_%u_phi_%.3f_repeat_%d_results.mat', P, phi, j);
        individual_filepath = fullfile(individual_data_dir, individual_filename);

        if isfile(individual_filepath)
            loaded = load(individual_filepath, 'run_data');
            if isfield(loaded.run_data, 'residual') && ~isempty(loaded.run_data.residual)
                residuals_all(i, j) = loaded.run_data.residual;
            elseif isfield(loaded.run_data, 'sol')
                % Fallback: extract from sol struct
                sol = loaded.run_data.sol;
                if isfield(sol, 'body_rel_res_max') && ~isempty(sol.body_rel_res_max)
                    residuals_all(i, j) = max(sol.body_rel_res_max);
                elseif isfield(sol, 'rel_res')
                    residuals_all(i, j) = sol.rel_res;
                end
            end
        else
            warning('Could not find file: %s', individual_filepath);
        end
    end
end
end

%% Helper function to extract n_close from individual files
function n_close_all = extract_n_close_from_files(data_root, P_range, n_repeats)
% Extract n_close from individual saved run data files
% Used as fallback when n_close was not stored in aggregated results
n_p = length(P_range);
n_close_all = zeros(n_p, n_repeats);
individual_data_dir = fullfile(data_root, 'random_close_packing');
phi = 0.65;  % From script parameters

for i = 1:n_p
    P = P_range(i);
    for j = 1:n_repeats
        individual_filename = sprintf('P_%u_phi_%.3f_repeat_%d_results.mat', P, phi, j);
        individual_filepath = fullfile(individual_data_dir, individual_filename);

        if isfile(individual_filepath)
            loaded = load(individual_filepath, 'run_data');
            if isfield(loaded.run_data, 'n_close')
                n_close_all(i, j) = loaded.run_data.n_close;
            else
                % Recompute if not available
                q = loaded.run_data.q;
                delta_pair = 0.2;
                rad = 1;
                n_close_all(i, j) = count_close_pairs(q, delta_pair, rad);
            end
        else
            warning('Could not find file: %s', individual_filepath);
        end
    end
end
end

%% Utility function for filled area in linear plot
function fill_between_wrapper(ax, x, y_min, y_max, alpha)
% Helper to fill between two curves (linear scale)
x_fill = [x(:)', fliplr(x(:)')];
y_fill = [y_min(:)', fliplr(y_max(:)')];
fill(ax, x_fill, y_fill, 'b', 'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

%% Utility function for filled area in log plot
function fill_between_log_wrapper(ax, x, y_min, y_max, alpha)
% Helper to fill between two curves (log scale)
x_fill = [x(:)', fliplr(x(:)')];
y_fill = [y_min(:)', fliplr(y_max(:)')];
fill(ax, x_fill, y_fill, 'b', 'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
