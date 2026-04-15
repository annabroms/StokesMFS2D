% Study how RAM usage scales with particle count for the peanut-compressed
% Stokes mobility solver on random Monte Carlo disc packings at phi = 0.55.

close all;
clear;

script_name = mfilename;
script_date = 'Apr 15, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root,'startup.m'));
end

fprintf('=== %s (%s) ===\n', script_name, script_date);
fprintf('Measurement mode: fresh external MATLAB process per (P, run).\n');
fprintf(['Default range is capped conservatively for desktop MATLAB, where ', ...
    'the parent session and child batch job coexist in RAM.\n']);

phi = 0.55;

if ~exist('P_values','var') || isempty(P_values)
    P_values = [100; 200; 400]; 
    %P_values = 200;
end

if ~exist('n_runs','var') || isempty(n_runs)
    n_runs = 3;
end

if ~exist('results_path','var') || isempty(results_path)
    results_path = fullfile(repo_root,'experiments','apr15_RAM_dependence_results.mat');
end

nP = numel(P_values);
baseline_bytes = nan(nP,n_runs);
overall_peak_bytes = nan(nP,n_runs);
overall_delta_bytes = nan(nP,n_runs);
precomp_peak_bytes = nan(nP,n_runs);
solve_peak_bytes = nan(nP,n_runs);
postprocess_peak_bytes = nan(nP,n_runs);
precomp_delta_bytes = nan(nP,n_runs);
solve_delta_bytes = nan(nP,n_runs);
postprocess_delta_bytes = nan(nP,n_runs);
solve_times = nan(nP,n_runs);
gmres_iters = nan(nP,n_runs);
gmres_unknowns = nan(nP,n_runs);
rel_residuals = nan(nP,n_runs);

run_data = repmat(struct( ...
    'P',nan, ...
    'run_index',nan, ...
    'geom_seed',nan, ...
    'geom_meta',struct(), ...
    'it',nan, ...
    'gmres_unknowns',nan, ...
    'rel_res',nan, ...
    'solve_time',nan, ...
    'baseline_bytes',nan, ...
    'overall_peak_bytes',nan, ...
    'overall_delta_bytes',nan, ...
    'precomp_peak_bytes',nan, ...
    'solve_peak_bytes',nan, ...
    'postprocess_peak_bytes',nan, ...
    'precomp_delta_bytes',nan, ...
    'solve_delta_bytes',nan, ...
    'postprocess_delta_bytes',nan), nP,n_runs);

save_progress(results_path,script_name,script_date,phi,P_values,n_runs, ...
    baseline_bytes,overall_peak_bytes,overall_delta_bytes, ...
    precomp_peak_bytes,solve_peak_bytes,postprocess_peak_bytes, ...
    precomp_delta_bytes,solve_delta_bytes,postprocess_delta_bytes, ...
    solve_times,gmres_iters,gmres_unknowns,rel_residuals,run_data);

for ip = 1:nP
    P = P_values(ip);
    fprintf('\nP = %d\n', P);

    for irun = 1:n_runs
        try
            result = run_case_external(repo_root,P,irun);
        catch ME
            save_progress(results_path,script_name,script_date,phi,P_values,n_runs, ...
                baseline_bytes,overall_peak_bytes,overall_delta_bytes, ...
                precomp_peak_bytes,solve_peak_bytes,postprocess_peak_bytes, ...
                precomp_delta_bytes,solve_delta_bytes,postprocess_delta_bytes, ...
                solve_times,gmres_iters,gmres_unknowns,rel_residuals,run_data);
            fprintf(2,'Partial results were saved to:\n  %s\n',results_path);
            rethrow(ME)
        end

        baseline_bytes(ip,irun) = result.baseline_bytes;
        overall_peak_bytes(ip,irun) = result.overall_peak_bytes;
        overall_delta_bytes(ip,irun) = result.overall_delta_bytes;
        precomp_peak_bytes(ip,irun) = result.precomp_peak_bytes;
        solve_peak_bytes(ip,irun) = result.solve_peak_bytes;
        postprocess_peak_bytes(ip,irun) = result.postprocess_peak_bytes;
        precomp_delta_bytes(ip,irun) = result.precomp_delta_bytes;
        solve_delta_bytes(ip,irun) = result.solve_delta_bytes;
        postprocess_delta_bytes(ip,irun) = result.postprocess_delta_bytes;
        solve_times(ip,irun) = result.solve_time;
        gmres_iters(ip,irun) = result.it;
        gmres_unknowns(ip,irun) = result.gmres_unknowns;
        rel_residuals(ip,irun) = result.rel_res;
        run_data(ip,irun).P = result.P;
        run_data(ip,irun).run_index = result.run_index;
        run_data(ip,irun).geom_seed = result.geom_seed;
        run_data(ip,irun).geom_meta = result.geom_meta;
        run_data(ip,irun).it = result.it;
        run_data(ip,irun).gmres_unknowns = result.gmres_unknowns;
        run_data(ip,irun).rel_res = result.rel_res;
        run_data(ip,irun).solve_time = result.solve_time;
        run_data(ip,irun).baseline_bytes = result.baseline_bytes;
        run_data(ip,irun).overall_peak_bytes = result.overall_peak_bytes;
        run_data(ip,irun).overall_delta_bytes = result.overall_delta_bytes;
        run_data(ip,irun).precomp_peak_bytes = result.precomp_peak_bytes;
        run_data(ip,irun).solve_peak_bytes = result.solve_peak_bytes;
        run_data(ip,irun).postprocess_peak_bytes = result.postprocess_peak_bytes;
        run_data(ip,irun).precomp_delta_bytes = result.precomp_delta_bytes;
        run_data(ip,irun).solve_delta_bytes = result.solve_delta_bytes;
        run_data(ip,irun).postprocess_delta_bytes = result.postprocess_delta_bytes;

        fprintf(['  run %d/%d: it=%d, unknowns=%d, rel_res=%.3e, time=%.2fs, ', ...
            'baseline=%.3f GiB, peak=%.3f GiB, delta=%.3f GiB\n'], ...
            irun,n_runs,result.it,result.gmres_unknowns,result.rel_res, ...
            result.solve_time,to_gib(result.baseline_bytes), ...
            to_gib(result.overall_peak_bytes),to_gib(result.overall_delta_bytes));

        save_progress(results_path,script_name,script_date,phi,P_values,n_runs, ...
            baseline_bytes,overall_peak_bytes,overall_delta_bytes, ...
            precomp_peak_bytes,solve_peak_bytes,postprocess_peak_bytes, ...
            precomp_delta_bytes,solve_delta_bytes,postprocess_delta_bytes, ...
            solve_times,gmres_iters,gmres_unknowns,rel_residuals,run_data);
    end
end

baseline_gib = to_gib(baseline_bytes);
overall_peak_gib = to_gib(overall_peak_bytes);
overall_delta_gib = to_gib(overall_delta_bytes);
precomp_peak_gib = to_gib(precomp_peak_bytes);
solve_peak_gib = to_gib(solve_peak_bytes);
postprocess_peak_gib = to_gib(postprocess_peak_bytes);
precomp_delta_gib = to_gib(precomp_delta_bytes);
solve_delta_gib = to_gib(solve_delta_bytes);
postprocess_delta_gib = to_gib(postprocess_delta_bytes);

fprintf('\nSummary (GiB, mean over runs):\n');
fprintf(['  %-5s %-12s %-12s %-12s %-12s %-12s %-12s\n'], ...
    'P','baseline','peak_all','delta_all','peak_pre','peak_sol','peak_post');
for ip = 1:nP
    fprintf('  %-5d %-12.3f %-12.3f %-12.3f %-12.3f %-12.3f %-12.3f\n', ...
        P_values(ip), ...
        mean(baseline_gib(ip,:),2,'omitnan'), ...
        mean(overall_peak_gib(ip,:),2,'omitnan'), ...
        mean(overall_delta_gib(ip,:),2,'omitnan'), ...
        mean(precomp_peak_gib(ip,:),2,'omitnan'), ...
        mean(solve_peak_gib(ip,:),2,'omitnan'), ...
        mean(postprocess_peak_gib(ip,:),2,'omitnan'));
end

colors = struct();
colors.overall = [0.10 0.10 0.10];
colors.precomp = [0.00 0.45 0.74];
colors.solve = [0.85 0.33 0.10];
colors.postprocess = [0.47 0.67 0.19];

figure('Name','apr15 RAM dependence','Color','w');
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
hold on;
plot_metric_series(P_values,precomp_peak_gib,'precomp peak',colors.precomp,n_runs);
plot_metric_series(P_values,solve_peak_gib,'solve peak',colors.solve,n_runs);
plot_metric_series(P_values,postprocess_peak_gib,'postprocess peak',colors.postprocess,n_runs);
plot_metric_series(P_values,overall_peak_gib,'overall peak',colors.overall,n_runs);
grid on;
xlabel('P');
ylabel('Peak MATLAB process RAM (GiB)');
title(sprintf('Peak MATLAB process RAM vs P at \\phi = %.2f', phi), ...
    'Interpreter','none');
legend('Location','northwest');

nexttile;
hold on;
plot_metric_series(P_values,precomp_delta_gib,'precomp delta',colors.precomp,n_runs);
plot_metric_series(P_values,solve_delta_gib,'solve delta',colors.solve,n_runs);
plot_metric_series(P_values,postprocess_delta_gib,'postprocess delta',colors.postprocess,n_runs);
plot_metric_series(P_values,overall_delta_gib,'overall delta',colors.overall,n_runs);
grid on;
xlabel('P');
ylabel('RAM increase above solver-start baseline (GiB)');
title('RAM increase above solver-start baseline vs P','Interpreter','none');
legend('Location','northwest');

title(t,'Peanut mobility RAM dependence','Interpreter','none');

save_progress(results_path,script_name,script_date,phi,P_values,n_runs, ...
    baseline_bytes,overall_peak_bytes,overall_delta_bytes, ...
    precomp_peak_bytes,solve_peak_bytes,postprocess_peak_bytes, ...
    precomp_delta_bytes,solve_delta_bytes,postprocess_delta_bytes, ...
    solve_times,gmres_iters,gmres_unknowns,rel_residuals,run_data);
fprintf('Saved results to:\n  %s\n',results_path);

function result = run_case_external(repo_root,P,run_index)
output_path = [tempname(tempdir), '.mat'];
cleanup_obj = onCleanup(@() delete_if_exists(output_path)); %#ok<NASGU>

batch_expr = sprintf([ ...
    'cd(''%s''); ' ...
    'run(fullfile(''%s'',''startup.m'')); ' ...
    'set(0,''DefaultFigureVisible'',''off''); ' ...
    'apr15_RAM_dependence_case(%d,%d,''%s'');'], ...
    escape_matlab_string(repo_root), ...
    escape_matlab_string(repo_root), ...
    P, ...
    run_index, ...
    escape_matlab_string(output_path));

cmd = sprintf([ ...
    'env OMP_NUM_THREADS=1 KMP_INIT_AT_FORK=FALSE matlab -batch "%s"'], ...
    escape_shell_double_quotes(batch_expr));

fprintf('  launching external MATLAB for P=%d, run=%d...\n',P,run_index);
[status,cmdout] = system(cmd);
if ~isempty(strtrim(cmdout))
    fprintf('%s',cmdout);
    if cmdout(end) ~= newline
        fprintf('\n');
    end
end

if status ~= 0
    error('apr15_RAM_dependence:externalRunFailed', ...
        'External MATLAB run failed for P=%d, run=%d.',P,run_index);
end

loaded = load(output_path,'result');
result = loaded.result;
end

function text = escape_matlab_string(text)
text = strrep(text,'''','''''');
end

function text = escape_shell_double_quotes(text)
text = strrep(text,'\','\\');
text = strrep(text,'"','\"');
text = strrep(text,'$','\$');
text = strrep(text,'`','\`');
end

function delete_if_exists(path_to_file)
if exist(path_to_file,'file') == 2
    delete(path_to_file);
end
end

function save_progress(results_path,script_name,script_date,phi,P_values,n_runs, ...
        baseline_bytes,overall_peak_bytes,overall_delta_bytes, ...
        precomp_peak_bytes,solve_peak_bytes,postprocess_peak_bytes, ...
        precomp_delta_bytes,solve_delta_bytes,postprocess_delta_bytes, ...
        solve_times,gmres_iters,gmres_unknowns,rel_residuals,run_data)
save(results_path, ...
    'script_name','script_date','phi','P_values','n_runs', ...
    'baseline_bytes','overall_peak_bytes','overall_delta_bytes', ...
    'precomp_peak_bytes','solve_peak_bytes','postprocess_peak_bytes', ...
    'precomp_delta_bytes','solve_delta_bytes','postprocess_delta_bytes', ...
    'solve_times','gmres_iters','gmres_unknowns','rel_residuals','run_data');
end

function plot_metric_series(x,values,label,color,n_runs)
if n_runs == 1
    plot(x,values(:),'-o', ...
        'Color',color, ...
        'LineWidth',1.5, ...
        'MarkerFaceColor',color, ...
        'DisplayName',label);
    return
end

vmin = min(values,[],2,'omitnan');
vmax = max(values,[],2,'omitnan');
vmean = mean(values,2,'omitnan');
plot_band(x,vmin,vmax,color);
plot(x,vmean,'-o', ...
    'Color',color, ...
    'LineWidth',1.5, ...
    'MarkerFaceColor',color, ...
    'DisplayName',label);
end

function plot_band(x,ymin,ymax,color)
xx = [x(:); flipud(x(:))];
yy = [ymin(:); flipud(ymax(:))];
patch(xx,yy,color, ...
    'FaceAlpha',0.15, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');
end

function values_gib = to_gib(values_bytes)
values_gib = values_bytes/(1024^3);
end
