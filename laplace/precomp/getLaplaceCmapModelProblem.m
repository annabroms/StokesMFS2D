function problem = getLaplaceCmapModelProblem(model)
%GETLAPLACECMAPMODELPROBLEM Return the problem represented by a Cmap model.

if ~isstruct(model)
    error('getLaplaceCmapModelProblem:BadModel', ...
        'The interpolation model must be a struct.');
end
kind = getOptField(model,'kind','');
switch kind
    case 'laplace_capacitance_cmap_alpha'
        kind_problem = 'capacitance';
    case 'laplace_elastance_cmap_alpha'
        kind_problem = 'elastance';
    otherwise
        error('getLaplaceCmapModelProblem:BadModelKind', ...
            'Unknown Laplace Cmap interpolation model kind.');
end
if isfield(model,'problem') && ~isempty(model.problem)
    problem = resolveLaplaceInterpolationProblem( ...
        model.problem,mfilename);
    if ~strcmp(problem,kind_problem)
        error('getLaplaceCmapModelProblem:InconsistentProblem', ...
            'The model problem and model kind are inconsistent.');
    end
else
    problem = kind_problem;
end
end
