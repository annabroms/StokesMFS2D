function [state,cleanup_guard] = startRamCheck(opt,solver_name)
%STARTRAMCHECK Start optional memorygraph-based RAM tracking for a solver.
%
% Syntax:
%   [state,cleanup_guard] = startRamCheck(opt,solver_name)
%
% Inputs:
%   opt         - Solver options struct. Uses opt.RAM_check when present.
%   solver_name - Name shown in the RAM summary output.
%
% Outputs:
%   state         - RAM tracking state passed to markRamCheckPhase and
%                   finishRamCheck.
%   cleanup_guard - onCleanup handle that stops memorygraph on exit.
%
% Notes:
%   memorygraph samples MATLAB RAM asynchronously via `top`, so we pause
%   briefly after starting it before the solver makes large allocations.
%   The actual baseline is taken from the first recorded sample in the
%   final trace.
%
% See also: markRamCheckPhase, finishRamCheck, memorygraph.

sampling_dt = 0.1;

state = init_ram_state();
state.enabled = logical(getOptField(opt,'RAM_check',false));
state.sampling_dt = sampling_dt;
state.solver_name = solver_name;
cleanup_guard = [];

if ~state.enabled
    return
end

opts = struct();
opts.dt = sampling_dt;
memorygraph('start',opts);
cleanup_guard = onCleanup(@stop_memorygraph_safely);

% Let memorygraph record the initial MATLAB footprint before the solver
% enters the precomputation phase.
pause(max(5*sampling_dt,0.5));
end

function state = init_ram_state()
state = struct();
state.enabled = false;
state.sampling_dt = nan;
state.solver_name = '';
state.baseline_sample = 1;
state.baseline_bytes = nan;
state.phase_marks = struct('start',1);
end

function stop_memorygraph_safely()
try
    memorygraph('done');
catch
end
end
