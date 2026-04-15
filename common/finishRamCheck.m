function summary = finishRamCheck(state)
%FINISHRAMCHECK Finalize a RAM trace and summarize solver phases.
%
% Syntax:
%   summary = finishRamCheck(state)
%
% Input:
%   state - State returned by startRamCheck and updated by
%           markRamCheckPhase.
%
% Output:
%   summary - Struct with raw-byte RAM estimates for the full solve and
%             for the precompute, solve, and postprocess phases.
%
% See also: startRamCheck, markRamCheckPhase.

sampling_dt = get_sampling_dt(state);
summary = init_ram_summary(false,sampling_dt);

if ~is_ram_check_enabled(state)
    return
end

bytes = get_trace_bytes();
summary = init_ram_summary(true,state.sampling_dt);

if isempty(bytes)
    return
end

baseline_bytes = get_baseline_bytes(state,bytes);
n_samples = numel(bytes);
phase_bounds = get_phase_bounds(state,n_samples);

summary.baseline_bytes = baseline_bytes;
summary.overall_peak_bytes = max(bytes);
summary.overall_delta_bytes = max(bytes - baseline_bytes);

summary.precomp = build_phase_summary(bytes,baseline_bytes, ...
    phase_bounds.precomp_start,phase_bounds.precomp_end);
summary.solve = build_phase_summary(bytes,baseline_bytes, ...
    phase_bounds.solve_start,phase_bounds.solve_end);
summary.postprocess = build_phase_summary(bytes,baseline_bytes, ...
    phase_bounds.postprocess_start,phase_bounds.postprocess_end);

print_ram_summary(state.solver_name,summary);
end

function tf = is_ram_check_enabled(state)
tf = isstruct(state) && isfield(state,'enabled') && state.enabled;
end

function value = get_sampling_dt(state)
if isstruct(state) && isfield(state,'sampling_dt') && isfinite(state.sampling_dt)
    value = state.sampling_dt;
else
    value = 0.1;
end
end

function summary = init_ram_summary(enabled,sampling_dt)
phase = struct('peak_bytes',nan,'delta_bytes',nan, ...
    'start_sample',nan,'end_sample',nan);
summary = struct();
summary.enabled = enabled;
summary.sampling_dt = sampling_dt;
summary.baseline_bytes = nan;
summary.overall_peak_bytes = nan;
summary.overall_delta_bytes = nan;
summary.precomp = phase;
summary.solve = phase;
summary.postprocess = phase;
end

function bytes = get_trace_bytes()
[bytes,~,~,~,~,~] = memorygraph('get');
end

function baseline_bytes = get_baseline_bytes(state,bytes)
baseline_bytes = state.baseline_bytes;

if ~isfinite(baseline_bytes)
    baseline_bytes = bytes(1);
end
end

function bounds = get_phase_bounds(state,n_samples)
precomp_start = get_phase_mark(state,'start',1);
precomp_end = get_phase_mark(state,'precomp_end',n_samples);
solve_end = get_phase_mark(state,'solve_end',n_samples);

bounds = struct();
bounds.precomp_start = precomp_start;
bounds.precomp_end = precomp_end;
bounds.solve_start = precomp_end;
bounds.solve_end = solve_end;
bounds.postprocess_start = solve_end;
bounds.postprocess_end = n_samples;
end

function value = get_phase_mark(state,field_name,default_value)
value = default_value;

if isstruct(state) && isfield(state,'phase_marks') && ...
        isstruct(state.phase_marks) && isfield(state.phase_marks,field_name)
    value = state.phase_marks.(field_name);
elseif strcmp(field_name,'start') && isstruct(state) && isfield(state,'baseline_sample')
    value = state.baseline_sample;
end
end

function phase = build_phase_summary(bytes,baseline_bytes,start_sample,end_sample)
phase = struct('peak_bytes',nan,'delta_bytes',nan, ...
    'start_sample',nan,'end_sample',nan);

if isempty(bytes)
    return
end

n_samples = numel(bytes);
start_sample = max(1,min(n_samples,round(start_sample)));
end_sample = max(start_sample,min(n_samples,round(end_sample)));
segment = bytes(start_sample:end_sample);
peak_bytes = max(segment);

phase.peak_bytes = peak_bytes;
phase.delta_bytes = max(peak_bytes - baseline_bytes,0);
phase.start_sample = start_sample;
phase.end_sample = end_sample;
end

function print_ram_summary(solver_name,summary)
fprintf(['RAM estimate (%s): baseline %.3f GiB, overall peak %.3f GiB, ', ...
    'delta %.3f GiB\n'], ...
    solver_name, ...
    to_gib(summary.baseline_bytes), ...
    to_gib(summary.overall_peak_bytes), ...
    to_gib(summary.overall_delta_bytes));
fprintf('  precomp     peak %.3f GiB, delta %.3f GiB\n', ...
    to_gib(summary.precomp.peak_bytes),to_gib(summary.precomp.delta_bytes));
fprintf('  solve       peak %.3f GiB, delta %.3f GiB\n', ...
    to_gib(summary.solve.peak_bytes),to_gib(summary.solve.delta_bytes));
fprintf('  postprocess peak %.3f GiB, delta %.3f GiB\n', ...
    to_gib(summary.postprocess.peak_bytes),to_gib(summary.postprocess.delta_bytes));
end

function value = to_gib(bytes)
value = bytes/(1024^3);
end
