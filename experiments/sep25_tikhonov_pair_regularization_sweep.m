%% Tikhonov parameter sweep for peanut-compressed two-body corrections
% Compare the legacy truncated-SVD (TSVD) pair correction with zero-order
% Tikhonov regularization for four problems on one reproducible random
% close packing of ten disks:
%
%   1. Laplace capacitance: prescribed body voltages, recovered charges.
%   2. Laplace elastance: prescribed zero-total body charges, voltages.
%   3. Stokes resistance: prescribed rigid velocities, forces/torques.
%   4. Stokes mobility: prescribed zero-total force and random torques,
%      recovered rigid velocities.
%
% The Stokes geometry has unit-radius disks. The Laplace geometry and box
% are scaled by laplace_radius=2, avoiding the unit logarithmic-capacity
% special case while retaining exactly the same dimensionless packing.
%
% For use_tikhonov=true, tikhonov_tol is the dimensionless parameter
%
%                  eta = lambda/sigma_max(A),
%
% and each pair or peanut pseudoinverse applies the spectral filter
%
%                  sigma/(sigma^2 + lambda^2).
%
% The one-body preconditioners remain unchanged. Low-rank Cmap compression
% is disabled, so the comparison isolates regularization of the two-body
% pseudoinverses. The script reports each solution's change from TSVD,
% independent-boundary residual, GMRES iterations, and wall time. The
% largest eta satisfying the configurable common acceptance criteria is
% printed as a starting recommendation; the full tables and plots should
% be used to choose a tolerance appropriate to the desired output accuracy.
%
% Anna Broms, Sep 2026

close all;

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s ===\n',mfilename);

%% Sweep and acceptance settings
% The first four values resolve the neighborhood of the legacy Laplace
% cutoff 1e-14. The final two probe the larger legacy Stokes fine-pair
% cutoff (1e-11). Add or remove values here for a narrower follow-up run.
if ~exist('tikhonov_tol_values','var') || isempty(tikhonov_tol_values)
    tikhonov_tol_values = [3e-15 1e-14 3e-14 1e-13 1e-12 1e-11];
end
if ~exist('output_change_tol','var') || isempty(output_change_tol)
    output_change_tol = 1e-4;
end
if ~exist('residual_growth_limit','var') || isempty(residual_growth_limit)
    residual_growth_limit = 10;
end
if ~exist('plotfig','var') || isempty(plotfig)
    plotfig = true;
end
if ~exist('save_results','var') || isempty(save_results)
    save_results = false;
end

%% Reproducible random close packing and discretisation
P = 10;
stokes_radius = 1;
laplace_radius = 2;
packing_fraction = 0.65;
minimum_gap = 1e-3*stokes_radius;
pair_gap = 0.2*stokes_radius;
geometry_seed = 250925;
load_seed = 250926;

N_c = 60;
N_f = 60;
N_peanut = 240;
Nclust = 80;
gmres_tol = 1e-8;
maxit = 500;

geom_opt = struct('domain','boxed','phi',packing_fraction, ...
    'rad',stokes_radius,'min_gap',minimum_gap,'n_sweeps',200, ...
    'rng_seed',geometry_seed,'visualise',false);
[q_stokes,geometry_meta] = random_discs_mc(P,geom_opt);
q_laplace = laplace_radius*q_stokes/stokes_radius;

[n_pairs,pairs,gaps] = count_close_pairs(q_stokes,pair_gap,stokes_radius);
if n_pairs == 0
    error('sep25_tikhonov_pair_regularization_sweep:NoClosePairs', ...
        'The generated packing has no pairs below the correction threshold.');
end

fprintf(['Geometry: P=%d, phi=%.3f, measured min gap=%.3e, ', ...
    'close pairs=%d below %.2f\n'], ...
    P,geometry_meta.phi,geometry_meta.min_surface_gap,n_pairs,pair_gap);
fprintf('Close-pair gap range: [%.3e, %.3e]\n',min(gaps),max(gaps));
fprintf('N_c=%d, N_f=%d, N_peanut=%d, Nclust=%d, gmres_tol=%.1e\n\n', ...
    N_c,N_f,N_peanut,Nclust,gmres_tol);

%% Fixed random loads for all regularization parameters
rng(load_seed,'twister');
loads = struct();
loads.voltage = randn(P,1);
loads.charge = randn(P,1);
loads.charge = loads.charge-mean(loads.charge); % exterior elastance compatibility
loads.velocity = randn(P,2);
loads.angular_velocity = randn(P,1);
loads.force = randn(P,2);
loads.force = loads.force-mean(loads.force,1); % 2D mobility compatibility
loads.torque = randn(P,1);

laplace_opt = make_laplace_options(P,laplace_radius,N_c,N_f, ...
    N_peanut,Nclust,gmres_tol,pair_gap*laplace_radius/stokes_radius);
stokes_opt = make_stokes_options(P,stokes_radius,N_c,N_f, ...
    N_peanut,Nclust,gmres_tol,maxit,pair_gap);

%% Legacy TSVD reference
fprintf('--- Legacy TSVD reference ---\n');
laplace_opt.use_tikhonov = false;
stokes_opt.use_tikhonov = false;
baseline = run_four_solvers(q_laplace,q_stokes,loads,laplace_opt,stokes_opt);

%% Tikhonov sweep
n_tol = numel(tikhonov_tol_values);
runs = repmat(init_suite_result(),n_tol,1);
for k = 1:n_tol
    eta = tikhonov_tol_values(k);
    fprintf('\n--- Tikhonov %d/%d: eta=lambda/sigma_max=%.3e ---\n', ...
        k,n_tol,eta);
    laplace_run_opt = laplace_opt;
    stokes_run_opt = stokes_opt;
    laplace_run_opt.use_tikhonov = true;
    stokes_run_opt.use_tikhonov = true;
    laplace_run_opt.tikhonov_tol = eta;
    stokes_run_opt.tikhonov_tol = eta;
    runs(k) = run_four_solvers(q_laplace,q_stokes,loads, ...
        laplace_run_opt,stokes_run_opt);
    runs(k).tikhonov_tol = eta;
    runs(k) = add_reference_changes(runs(k),baseline);
end

%% Report tables and a conservative common recommendation
problem_names = {'capacitance','elastance','resistance','mobility'};
fprintf('\nRelative output change from the legacy TSVD solve:\n');
fprintf('%12s %14s %14s %14s %14s\n', ...
    'eta','capacitance','elastance','resistance','mobility');
for k = 1:n_tol
    fprintf('%12.3e %14.3e %14.3e %14.3e %14.3e\n', ...
        runs(k).tikhonov_tol,runs(k).capacitance.output_change, ...
        runs(k).elastance.output_change,runs(k).resistance.output_change, ...
        runs(k).mobility.output_change);
end

fprintf('\nIndependent-boundary residuals:\n');
fprintf('%12s %14s %14s %14s %14s\n', ...
    'eta','capacitance','elastance','resistance','mobility');
fprintf('%12s %14.3e %14.3e %14.3e %14.3e\n', ...
    'TSVD',baseline.capacitance.residual,baseline.elastance.residual, ...
    baseline.resistance.residual,baseline.mobility.residual);
for k = 1:n_tol
    fprintf('%12.3e %14.3e %14.3e %14.3e %14.3e\n', ...
        runs(k).tikhonov_tol,runs(k).capacitance.residual, ...
        runs(k).elastance.residual,runs(k).resistance.residual, ...
        runs(k).mobility.residual);
end

fprintf('\nGMRES iterations:\n');
fprintf('%12s %14s %14s %14s %14s\n', ...
    'eta','capacitance','elastance','resistance','mobility');
fprintf('%12s %14d %14d %14d %14d\n','TSVD', ...
    baseline.capacitance.iterations,baseline.elastance.iterations, ...
    baseline.resistance.iterations,baseline.mobility.iterations);
for k = 1:n_tol
    fprintf('%12.3e %14d %14d %14d %14d\n', ...
        runs(k).tikhonov_tol,runs(k).capacitance.iterations, ...
        runs(k).elastance.iterations,runs(k).resistance.iterations, ...
        runs(k).mobility.iterations);
end

accepted = false(n_tol,numel(problem_names));
for k = 1:n_tol
    for j = 1:numel(problem_names)
        name = problem_names{j};
        ref_residual = baseline.(name).residual;
        residual_limit = max(residual_growth_limit*ref_residual,10*gmres_tol);
        accepted(k,j) = runs(k).(name).output_change <= output_change_tol && ...
            runs(k).(name).residual <= residual_limit;
    end
end
common_accepted = all(accepted,2);
recommended_index = find(common_accepted,1,'last');
fprintf(['\nAcceptance rule: output change <= %.1e and residual <= max(', ...
    '%.1f*TSVD residual, %.1e).\n'], ...
    output_change_tol,residual_growth_limit,10*gmres_tol);
if isempty(recommended_index)
    fprintf(['No swept eta satisfies the rule for all four problems. ', ...
        'Inspect the per-problem tables and reduce eta.\n']);
    recommended_tikhonov_tol = NaN;
else
    recommended_tikhonov_tol = tikhonov_tol_values(recommended_index);
    fprintf('Largest common accepted eta: %.3e\n',recommended_tikhonov_tol);
end

%% Plots
if plotfig
    labels = {'Laplace capacitance','Laplace elastance', ...
        'Stokes resistance','Stokes mobility'};
    colors = lines(4);

    figure('Name','Tikhonov output change','Color','w');
    hold on;
    for j = 1:numel(problem_names)
        vals = arrayfun(@(r) r.(problem_names{j}).output_change,runs);
        loglog(tikhonov_tol_values,vals,'o-','LineWidth',1.4, ...
            'Color',colors(j,:),'DisplayName',labels{j});
    end
    yline(output_change_tol,'k--','acceptance threshold', ...
        'HandleVisibility','off');
    xlabel('\eta = \lambda/\sigma_{max}');
    ylabel('relative output change from TSVD');
    title('Two-body Tikhonov regularization: solution sensitivity');
    legend('Location','best'); grid on;

    figure('Name','Tikhonov boundary residual','Color','w');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for j = 1:numel(problem_names)
        name = problem_names{j};
        nexttile;
        vals = arrayfun(@(r) r.(name).residual,runs);
        loglog(tikhonov_tol_values,vals,'o-','LineWidth',1.4);
        hold on;
        yline(baseline.(name).residual,'--','TSVD');
        xlabel('\eta = \lambda/\sigma_{max}');
        ylabel('independent-boundary residual');
        title(labels{j}); grid on;
    end

    figure('Name','Tikhonov GMRES iterations','Color','w');
    hold on;
    for j = 1:numel(problem_names)
        vals = arrayfun(@(r) r.(problem_names{j}).iterations,runs);
        semilogx(tikhonov_tol_values,vals,'o-','LineWidth',1.4, ...
            'Color',colors(j,:),'DisplayName',labels{j});
    end
    xlabel('\eta = \lambda/\sigma_{max}');
    ylabel('GMRES iterations');
    title('Solver convergence across the Tikhonov sweep');
    legend('Location','best'); grid on;

    plot_packing(q_stokes,stokes_radius,pairs,geometry_meta);
end

if save_results
    results_file = fullfile(fileparts(mfilename('fullpath')), ...
        'sep25_tikhonov_pair_regularization_sweep_results.mat');
    save(results_file,'baseline','runs','tikhonov_tol_values', ...
        'recommended_tikhonov_tol','accepted','q_stokes','q_laplace', ...
        'geometry_meta','loads','laplace_opt','stokes_opt', ...
        'output_change_tol','residual_growth_limit','gmres_tol');
    fprintf('Saved results to %s\n',results_file);
end

result = struct('baseline',baseline,'runs',runs, ...
    'tikhonov_tol_values',tikhonov_tol_values, ...
    'recommended_tikhonov_tol',recommended_tikhonov_tol, ...
    'accepted',accepted,'output_change_tol',output_change_tol, ...
    'residual_growth_limit',residual_growth_limit, ...
    'q_stokes',q_stokes,'q_laplace',q_laplace, ...
    'geometry_meta',geometry_meta,'loads',loads);

fprintf('\nSweep complete. Results are available in workspace variable result.\n');

function opt = make_laplace_options(P,R,N_c,N_f,N_peanut,Nclust, ...
    gmres_tol,delta_pair)
opt = getLaplace2Dparams(P,R,N_c,N_f);
opt.delta_pair = delta_pair;
opt.N_peanut = N_peanut;
opt.Nclust = Nclust;
opt.gmres_tol = gmres_tol;
opt.cmap = true;
opt.compress_cmap = false;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.use_big_sparse = false;
opt.use_fmm = false;
opt.get_bndry_field = true;
opt.get_precomp_time = true;
opt.get_solve_time = true;
opt.visualise_sol = false;
opt.visualise_grid = false;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.debug = false;
end

function opt = make_stokes_options(P,R,N_c,N_f,N_peanut,Nclust, ...
    gmres_tol,maxit,delta_pair)
opt = get2Dparams(P,N_c,N_f);
opt.rad = R;
opt.delta_pair = delta_pair;
opt.N_peanut = N_peanut;
opt.Nclust = Nclust;
opt.gmres_tol = gmres_tol;
opt.maxit = maxit;
opt.cmap = true;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.parallel_big_sparse_build = false;
opt.parallel_solve = false;
opt.use_big_sparse = false;
opt.use_fmm = false;
opt.self_correct = true;
opt.use_dense = true;
opt.get_bndry_field = true;
opt.get_precomp_time = true;
opt.get_solve_time = true;
opt.visualise_sol = false;
opt.visualise_grid = false;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.debug = false;
opt.column_weight = false;
opt.left_weight = false;
opt.solve_threads = 1;
end

function suite = run_four_solvers(q_laplace,q_stokes,loads,laplace_opt,stokes_opt)
suite = init_suite_result();

timer = tic;
[output,sol] = solve_cap_peanut(q_laplace,loads.voltage,laplace_opt);
suite.capacitance = summarize_solve(output,sol,sol.maxres,toc(timer));
print_run_summary('Laplace capacitance',suite.capacitance);

timer = tic;
[output,sol] = solve_elast_peanut(q_laplace,loads.charge,laplace_opt);
suite.elastance = summarize_solve(output,sol,sol.maxres,toc(timer));
print_run_summary('Laplace elastance',suite.elastance);

timer = tic;
[output,sol] = solve_res_peanut_enhanced(q_stokes,loads.velocity, ...
    loads.angular_velocity,stokes_opt);
suite.resistance = summarize_solve(output,sol,sol.rel_res,toc(timer));
print_run_summary('Stokes resistance',suite.resistance);

timer = tic;
[output,sol] = solve_mob_peanut_enhanced(q_stokes,loads.force, ...
    loads.torque,stokes_opt);
suite.mobility = summarize_solve(output,sol,sol.rel_res,toc(timer));
print_run_summary('Stokes mobility',suite.mobility);
end

function summary = summarize_solve(output,sol,residual,wall_time)
summary = init_problem_result();
summary.output = output;
summary.residual = residual;
summary.iterations = sol.it;
summary.gmres_final = sol.resvec(end);
summary.wall_time = wall_time;
if isfield(sol,'precomp_time') && isfield(sol.precomp_time,'total')
    summary.precomp_time = sol.precomp_time.total;
end
if isfield(sol,'solve_time') && isfield(sol.solve_time,'total')
    summary.solve_time = sol.solve_time.total;
end
end

function print_run_summary(label,summary)
fprintf('  %-22s it=%4d residual=%.3e wall=%7.2fs\n', ...
    label,summary.iterations,summary.residual,summary.wall_time);
end

function suite = add_reference_changes(suite,reference)
names = {'capacitance','elastance','resistance','mobility'};
for j = 1:numel(names)
    name = names{j};
    suite.(name).output_change = relerr( ...
        suite.(name).output,reference.(name).output);
end
end

function suite = init_suite_result()
empty_problem = init_problem_result();
suite = struct('tikhonov_tol',NaN,'capacitance',empty_problem, ...
    'elastance',empty_problem,'resistance',empty_problem, ...
    'mobility',empty_problem);
end

function result = init_problem_result()
result = struct('output',[],'output_change',NaN,'residual',NaN, ...
    'iterations',0,'gmres_final',NaN,'wall_time',NaN, ...
    'precomp_time',NaN,'solve_time',NaN);
end

function value = relerr(actual,reference)
value = norm(actual(:)-reference(:),inf)/max(1,norm(reference(:),inf));
end

function plot_packing(q,rad,pairs,geometry_meta)
figure('Name','Tikhonov sweep packing','Color','w');
theta = linspace(0,2*pi,200);
hold on;
for k = 1:numel(q)
    z = q(k)+rad*exp(1i*theta);
    plot(real(z),imag(z),'k-','LineWidth',0.9);
    text(real(q(k)),imag(q(k)),sprintf('%d',k), ...
        'HorizontalAlignment','center');
end
for row = 1:size(pairs,1)
    z = q(pairs(row,:));
    plot(real(z),imag(z),'-','Color',[0.85 0.33 0.10], ...
        'LineWidth',1.2);
end
half_width = geometry_meta.L/2;
plot([-half_width half_width half_width -half_width -half_width], ...
    [-half_width -half_width half_width half_width -half_width], ...
    'k--','LineWidth',0.8);
axis equal; grid on;
xlabel('x'); ylabel('y');
title(sprintf('Random close packing: P=%d, close pairs=%d', ...
    numel(q),size(pairs,1)));
end
