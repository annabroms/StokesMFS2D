% Compare aggregate RAM usage for streamed low-RAM Stokes pair precompute.
%
% Each case runs in a fresh MATLAB process. The parent process samples the
% child MATLAB process tree, so local parallel workers are included.

close all;

script_name = mfilename;
script_date = 'Apr 20, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root,'startup.m'));
end

fprintf('=== %s (%s) ===\n',script_name,script_date);
fprintf('Measurement mode: fresh external MATLAB process per variant.\n');
fprintf('RAM measurement: aggregate RSS of child MATLAB plus descendants.\n');

if ~exist('P_values','var') || isempty(P_values)
    P_values = 100;
end
if ~exist('n_runs','var') || isempty(n_runs)
    n_runs = 1;
end
if ~exist('pool_sizes','var') || isempty(pool_sizes)
    pool_sizes = [2 4];
end
if ~exist('results_path','var') || isempty(results_path)
    results_path = fullfile(repo_root,'experiments', ...
        'apr20_parallel_precomp_ram_results.mat');
end

variants = { ...
    'serial_full', ...
    'parallel_streamed_full', ...
    'parallel_streamed_auto_slim', ...
    'parallel_streamed_auto_slim_nodense'};

results = cell(0,1);

for ip = 1:numel(P_values)
    P = P_values(ip);
    for irun = 1:n_runs
        for ivar = 1:numel(variants)
            variant = variants{ivar};
            if strcmp(variant,'serial_full')
                sizes_this = 0;
            else
                sizes_this = pool_sizes(:).';
            end

            for isize = 1:numel(sizes_this)
                pool_size = sizes_this(isize);
                fprintf('\nP=%d run=%d variant=%s pool=%d\n', ...
                    P,irun,variant,pool_size);
                result = run_case_external(repo_root,P,irun,variant,pool_size);
                results{end+1,1} = result; %#ok<SAGROW>
                print_result(result);
                save(results_path,'script_name','script_date','P_values', ...
                    'n_runs','pool_sizes','variants','results');
            end
        end
    end
end

fprintf('\n=== Summary ===\n');
fprintf(['%-36s %6s %8s %10s %10s %10s %10s %10s\n'], ...
    'variant','pool','payload','pre_GiB','peak_GiB','ret_GiB', ...
    'pair_s','total_s');
for k = 1:numel(results)
    print_summary_row(results{k});
end

save(results_path,'script_name','script_date','P_values','n_runs', ...
    'pool_sizes','variants','results');
fprintf('Saved results to:\n  %s\n',results_path);

function result = run_case_external(repo_root,P,run_index,variant,pool_size)
output_path = [tempname(tempdir), '.mat'];
monitor_path = [tempname(tempdir), '.txt'];
script_path = [tempname(tempdir), '.sh'];
cleanup_obj = onCleanup(@() delete_paths({output_path,monitor_path,script_path})); %#ok<NASGU>

batch_expr = sprintf([ ...
    'cd(''%s''); ' ...
    'run(fullfile(''%s'',''startup.m'')); ' ...
    'set(0,''DefaultFigureVisible'',''off''); ' ...
    'apr20_parallel_precomp_ram_case(%d,%d,''%s'',%d,''%s'');'], ...
    escape_matlab_string(repo_root), ...
    escape_matlab_string(repo_root), ...
    P, ...
    run_index, ...
    escape_matlab_string(variant), ...
    pool_size, ...
    escape_matlab_string(output_path));

write_monitor_script(script_path,monitor_path,batch_expr);
cmd = sprintf('bash "%s"',script_path);
tic;
[status,cmdout] = system(cmd);
external_time = toc;
if ~isempty(strtrim(cmdout))
    fprintf('%s',cmdout);
    if cmdout(end) ~= newline
        fprintf('\n');
    end
end
if status ~= 0
    error('apr20_parallel_precomp_ram_benchmark:externalRunFailed', ...
        'External MATLAB run failed for variant=%s, pool=%d.',variant,pool_size);
end

loaded = load(output_path,'result');
result = loaded.result;
monitor = load_monitor_samples(monitor_path);
result = add_aggregate_ram_metrics(result,monitor);
result.external_time = external_time;
end

function write_monitor_script(script_path,monitor_path,batch_expr)
fid = fopen(script_path,'w');
if fid < 0
    error('Could not create monitor script: %s',script_path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'#!/usr/bin/env bash\n');
fprintf(fid,'set +e\n');
fprintf(fid,'monitor_file="%s"\n',escape_shell_double_quotes(monitor_path));
fprintf(fid,'sum_related_rss() {\n');
fprintf(fid,'  root="$1"\n');
fprintf(fid,'  pgid="$2"\n');
fprintf(fid,'  pids=""\n');
fprintf(fid,'  if ps -p "$root" >/dev/null 2>&1; then\n');
fprintf(fid,'    pids="$root"\n');
fprintf(fid,'    frontier="$root"\n');
fprintf(fid,'    while [ -n "$frontier" ]; do\n');
fprintf(fid,'      children="$(ps -e -o pid= -o ppid= | awk -v parents="$frontier" ''BEGIN{split(parents,a," "); for (i in a) p[a[i]]=1} p[$2]{print $1}'')"\n');
fprintf(fid,'      if [ -z "$children" ]; then break; fi\n');
fprintf(fid,'      pids="$pids $children"\n');
fprintf(fid,'      frontier="$children"\n');
fprintf(fid,'    done\n');
fprintf(fid,'  fi\n');
fprintf(fid,'  if [ -n "$pgid" ]; then\n');
fprintf(fid,'    group_pids="$(ps -e -o pid= -o pgid= | awk -v pg="$pgid" ''$2==pg{print $1}'')"\n');
fprintf(fid,'    pids="$pids $group_pids"\n');
fprintf(fid,'  fi\n');
fprintf(fid,'  pids="$(printf "%%s\\n" $pids | awk ''NF && !seen[$1]++{print $1}'')"\n');
fprintf(fid,'  if [ -z "$pids" ]; then echo 0; return; fi\n');
fprintf(fid,'  pid_csv="$(printf "%%s\\n" $pids | paste -sd, -)"\n');
fprintf(fid,'  ps -o rss= -p "$pid_csv" 2>/dev/null | awk ''{s+=$1} END{printf "%%.0f\\n", s*1024}''\n');
fprintf(fid,'}\n');
fprintf(fid,'setsid env OMP_NUM_THREADS=1 KMP_INIT_AT_FORK=FALSE matlab -batch "%s" &\n', ...
    escape_shell_double_quotes(batch_expr));
fprintf(fid,'child=$!\n');
fprintf(fid,'child_pgid="$(ps -o pgid= -p "$child" | tr -d " " | tr -d "\\n")"\n');
fprintf(fid,'(\n');
fprintf(fid,'  while kill -0 "$child" 2>/dev/null || [ -n "$(ps -e -o pgid= | awk -v pg="$child_pgid" ''$1==pg{print; exit}'')" ]; do\n');
fprintf(fid,'    bytes="$(sum_related_rss "$child" "$child_pgid")"\n');
fprintf(fid,'    ts="$(date +%%s.%%N)"\n');
fprintf(fid,'    printf "%%s %%s\\n" "$ts" "$bytes" >> "$monitor_file"\n');
fprintf(fid,'    sleep 0.1\n');
fprintf(fid,'  done\n');
fprintf(fid,') &\n');
fprintf(fid,'monitor_pid=$!\n');
fprintf(fid,'wait "$child"\n');
fprintf(fid,'status=$?\n');
fprintf(fid,'kill "$monitor_pid" 2>/dev/null\n');
fprintf(fid,'wait "$monitor_pid" 2>/dev/null\n');
fprintf(fid,'exit "$status"\n');
end

function monitor = load_monitor_samples(monitor_path)
if exist(monitor_path,'file') ~= 2
    monitor = zeros(0,2);
    return
end
monitor = readmatrix(monitor_path,'FileType','text');
if isempty(monitor)
    monitor = zeros(0,2);
elseif size(monitor,2) == 1
    monitor = [monitor zeros(size(monitor))];
end
end

function result = add_aggregate_ram_metrics(result,monitor)
if isempty(monitor)
    result.aggregate_baseline_bytes = NaN;
    result.aggregate_precomp_peak_bytes = NaN;
    result.aggregate_total_peak_bytes = NaN;
    result.aggregate_retained_after_precomp_bytes = NaN;
    result.aggregate_final_bytes = NaN;
    result.monitor_sample_count = 0;
    return
end

ts = monitor(:,1);
bytes = monitor(:,2);
valid = isfinite(ts) & isfinite(bytes) & bytes > 0;
ts = ts(valid);
bytes = bytes(valid);
if isempty(bytes)
    result.aggregate_baseline_bytes = NaN;
    result.aggregate_precomp_peak_bytes = NaN;
    result.aggregate_total_peak_bytes = NaN;
    result.aggregate_retained_after_precomp_bytes = NaN;
    result.aggregate_final_bytes = NaN;
    result.monitor_sample_count = 0;
    return
end
result.aggregate_baseline_bytes = bytes(1);
result.aggregate_final_bytes = bytes(end);
result.monitor_sample_count = numel(bytes);

precomp_total = get_precomp_total_time(result.precomp_time);
precomp_start = result.solve_start_posix;
precomp_end = precomp_start + precomp_total;
solve_end = result.solve_end_posix;

pre_mask = (ts >= precomp_start) & (ts <= precomp_end);
solve_mask = (ts >= precomp_start) & (ts <= solve_end);
if any(pre_mask)
    result.aggregate_precomp_peak_bytes = max(bytes(pre_mask));
else
    result.aggregate_precomp_peak_bytes = NaN;
end
if any(solve_mask)
    result.aggregate_total_peak_bytes = max(bytes(solve_mask));
else
    result.aggregate_total_peak_bytes = max(bytes);
end

after_pre = find(ts >= precomp_end,1,'first');
if isempty(after_pre)
    result.aggregate_retained_after_precomp_bytes = bytes(end);
else
    result.aggregate_retained_after_precomp_bytes = bytes(after_pre);
end

result.pair_basis_time = get_pair_basis_time(result.precomp_time);
result.precomp_total_time = precomp_total;
end

function value = get_pair_basis_time(precomp_time)
if isstruct(precomp_time) && isfield(precomp_time,'pair_basis')
    value = precomp_time.pair_basis;
else
    value = NaN;
end
end

function value = get_precomp_total_time(precomp_time)
if isstruct(precomp_time) && isfield(precomp_time,'total')
    value = precomp_time.total;
else
    value = NaN;
end
end

function print_result(result)
fprintf(['  payload=%s backend=%s max_inflight=%g pair_basis=%.2fs ', ...
    'total=%.2fs agg_pre=%.3f GiB agg_peak=%.3f GiB retained=%.3f GiB\n'], ...
    result.payload_mode,result.parallel_backend,result.max_inflight, ...
    result.pair_basis_time,result.solve_time, ...
    to_gib(result.aggregate_precomp_peak_bytes), ...
    to_gib(result.aggregate_total_peak_bytes), ...
    to_gib(result.aggregate_retained_after_precomp_bytes));
end

function print_summary_row(result)
fprintf('%-36s %6d %8s %10.3f %10.3f %10.3f %10.2f %10.2f\n', ...
    result.variant, ...
    result.pool_size_requested, ...
    result.payload_mode, ...
    to_gib(result.aggregate_precomp_peak_bytes), ...
    to_gib(result.aggregate_total_peak_bytes), ...
    to_gib(result.aggregate_retained_after_precomp_bytes), ...
    result.pair_basis_time, ...
    result.solve_time);
end

function text = escape_matlab_string(text)
text = strrep(char(text),'''','''''');
end

function text = escape_shell_double_quotes(text)
text = strrep(char(text),'\','\\');
text = strrep(text,'"','\"');
text = strrep(text,'$','\$');
text = strrep(text,'`','\`');
end

function delete_paths(paths)
for k = 1:numel(paths)
    if exist(paths{k},'file') == 2
        delete(paths{k});
    end
end
end

function value = to_gib(bytes)
value = bytes/(1024^3);
end
