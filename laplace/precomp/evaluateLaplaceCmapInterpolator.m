function payload = evaluateLaplaceCmapInterpolator(model,delta)
%EVALUATELAPLACECMAPINTERPOLATOR Evaluate one production alpha-panel model.

if ~isstruct(model) || ~isfield(model,'version') || model.version ~= 3 || ...
        ~isfield(model,'kind') || ~any(strcmp(model.kind,{ ...
        'laplace_capacitance_cmap_alpha', ...
        'laplace_elastance_cmap_alpha'}))
    error('evaluateLaplaceCmapInterpolator:BadModel', ...
        'Input is not a version-3 Laplace Cmap interpolator.');
end
problem = getLaplaceCmapModelProblem(model);
tolerance = 1e-12*max(1,model.delta_max);
if delta < model.delta_min-tolerance || delta > model.delta_max+tolerance
    error('evaluateLaplaceCmapInterpolator:GapOutOfRange', ...
        'Gap %.16g is outside [%.16g, %.16g].', ...
        delta,model.delta_min,model.delta_max);
end
delta = min(max(delta,model.delta_min),model.delta_max);
alpha = acosh(1+delta/(2*model.R));

panel_index = find(arrayfun(@(p) ...
    alpha >= p.alpha_lo-10*eps(max(1,abs(p.alpha_lo))) && ...
    alpha <= p.alpha_hi+10*eps(max(1,abs(p.alpha_hi))), ...
    model.panels),1,'first');
if isempty(panel_index)
    error('evaluateLaplaceCmapInterpolator:NoPanel', ...
        'No alpha panel contains gap %.16g.',delta);
end
panel = model.panels(panel_index);
coefficients = barycentricCoefficients(alpha,panel.alpha_nodes, ...
    panel.barycentric_weights);
charge_coefficients = barycentricCoefficients(alpha, ...
    panel.charge_alpha_nodes,panel.charge_barycentric_weights);

payload = struct('problem',problem,'mode',model.mode,'panel_index',panel_index, ...
    'coefficients',coefficients,'Cmap',[],'Cmap_QV',[], ...
    'charge_coefficients',charge_coefficients,'B',[]);
payload.Cmap_QV = combineSnapshots(panel.QV_nodes,charge_coefficients);
switch model.mode
    case 'full'
        payload.Cmap = combineSnapshots(panel.C_nodes,coefficients);
    case {'reduced','reduced_noconst'}
        payload.B = combineSnapshots(panel.B_nodes,coefficients);
    otherwise
        error('evaluateLaplaceCmapInterpolator:BadMode', ...
            'Unsupported model mode "%s".',model.mode);
end
end

function coefficients = barycentricCoefficients(x,nodes,weights)
difference = x-nodes;
[distance,index] = min(abs(difference));
if distance <= 1e-12*max(1,max(abs(nodes)))
    coefficients = zeros(numel(nodes),1);
    coefficients(index) = 1;
    return
end
terms = weights./difference;
coefficients = terms/sum(terms);
end

function value = combineSnapshots(snapshots,coefficients)
[n1,n2,q] = size(snapshots);
value = reshape(reshape(snapshots,n1*n2,q)*coefficients,n1,n2);
end
