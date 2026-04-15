function state = markRamCheckPhase(state,phase_name)
%MARKRAMCHECKPHASE Record a RAM sampling boundary for the current solver.
%
% Syntax:
%   state = markRamCheckPhase(state,phase_name)
%
% Inputs:
%   state      - State returned by startRamCheck.
%   phase_name - Boundary name, for example 'precomp_end' or 'solve_end'.
%
% Output:
%   state - Updated RAM tracking state.
%
% See also: startRamCheck, finishRamCheck.

if ~isstruct(state) || ~isfield(state,'enabled') || ~state.enabled
    return
end

memorygraph('label',phase_name);
[bytes,~,~,~,~,~] = memorygraph('get');
if isempty(bytes)
    state.phase_marks.(phase_name) = 1;
    return
end

state.phase_marks.(phase_name) = numel(bytes);
end
