function out = manageSolveTimeMeasurement(action,varargin)
%MANAGESOLVETIMEMEASUREMENT Track GMRES wall-clock time and FMM sub-time.
%
% Syntax:
%   token = manageSolveTimeMeasurement('start',enabled)
%   manageSolveTimeMeasurement('add_fmm',dt)
%   solve_time = manageSolveTimeMeasurement('finish',token)
%   manageSolveTimeMeasurement('reset')
%
% Notes:
%   The GMRES total time is measured by the solver around helsing_gmres.
%   The FMM time is accumulated only while a timing session is active.

persistent state
if isempty(state)
    state = init_state();
end

switch lower(action)
    case 'start'
        enabled = true;
        if nargin >= 2 && ~isempty(varargin{1})
            enabled = logical(varargin{1});
        end
        state = init_state();
        state.active = enabled;

        out = struct('enabled',enabled,'total_timer',[]);
        if enabled
            out.total_timer = tic;
        end

    case 'add_fmm'
        out = [];
        if nargin < 2 || isempty(varargin{1}) || ~state.active
            return
        end
        dt = varargin{1};
        if isnumeric(dt) && isscalar(dt) && isfinite(dt) && dt >= 0
            state.fmm_time = state.fmm_time + dt;
        end

    case 'finish'
        out = init_summary();
        if nargin < 2 || isempty(varargin{1}) || ~isstruct(varargin{1}) || ...
                ~isfield(varargin{1},'enabled') || ~varargin{1}.enabled
            state = init_state();
            return
        end

        token = varargin{1};
        out.total = toc(token.total_timer);
        out.fmm = state.fmm_time;
        state = init_state();

    case 'reset'
        out = [];
        state = init_state();

    otherwise
        error('manageSolveTimeMeasurement:UnknownAction', ...
            'Unknown action "%s".', action);
end

end

function state = init_state()
state = struct('active',false,'fmm_time',0);
end

function summary = init_summary()
summary = struct('total',nan,'fmm',nan);
end
