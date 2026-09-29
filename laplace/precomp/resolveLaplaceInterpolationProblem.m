function problem = resolveLaplaceInterpolationProblem(problem,caller_name)
%RESOLVELAPLACEINTERPOLATIONPROBLEM Validate the interpolation problem.

if nargin < 1 || isempty(problem)
    problem = 'capacitance';
end
if nargin < 2 || isempty(caller_name)
    caller_name = mfilename;
end
if isstring(problem)
    if ~isscalar(problem)
        error([caller_name ':BadProblem'], ...
            'The interpolation problem must be scalar text.');
    end
    problem = char(problem);
end
if ~ischar(problem)
    error([caller_name ':BadProblem'], ...
        'The interpolation problem must be ''capacitance'' or ''elastance''.');
end
problem = lower(strtrim(problem));
if ~any(strcmp(problem,{'capacitance','elastance'}))
    error([caller_name ':BadProblem'], ...
        'The interpolation problem must be ''capacitance'' or ''elastance''.');
end
end
