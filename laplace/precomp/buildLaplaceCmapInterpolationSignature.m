function signature = buildLaplaceCmapInterpolationSignature( ...
        opt,include_targets,problem)
%BUILDLAPLACECMAPINTERPOLATIONSIGNATURE Parameter identity for saved models.

if nargin < 2
    include_targets = true;
end
if nargin < 3
    problem = getOptField(opt,'interpolation_problem','capacitance');
end
problem = resolveLaplaceInterpolationProblem(problem,mfilename);
[opt,mode] = configureLaplaceCmapInterpolation(opt,problem);

signature = struct();
signature.format_version = 3;
signature.trainer_version = 1;
signature.mode = mode;
signature.problem = problem;
signature.project_charge = strcmp(problem,'elastance');
signature.rad = opt.rad;
signature.N_cmap = getOptField(opt,'N_cmap',opt.N_c);
signature.N_f = opt.N_f;
signature.N_peanut = opt.N_peanut;
signature.a_f = opt.a_f;
signature.Rp_c = opt.Rp_c;
signature.Rp_f = opt.Rp_f;
signature.delta_min = opt.smallest_delta;
signature.delta_max = opt.delta_pair;
signature.beta = opt.beta;
signature.Nclust = opt.Nclust;
signature.ellipse_constant = logical(opt.ellipse_constant);
signature.use_tikhonov = logical(opt.use_tikhonov);
signature.tikhonov_tol = opt.tikhonov_tol;
signature.compress_cmap = logical(opt.compress_cmap);
signature.coordinate = 'alpha';
signature.panel_count_candidates = getOptField(opt, ...
    'interpolation_panel_count_candidates',[1 2 4 8]);
signature.node_candidates = getOptField(opt, ...
    'interpolation_node_candidates',3:2:17);
signature.panel_count_candidates = ...
    signature.panel_count_candidates(:).';
signature.node_candidates = signature.node_candidates(:).';
signature.validation_nodes = getOptField(opt, ...
    'interpolation_validation_nodes',33);
signature.audit_nodes = 4;
signature.reference_location = 'delta_midpoint';
if include_targets
    signature.interpolation_tol = opt.interpolation_tol;
    signature.volt_charge_interp_tol = opt.volt_charge_interp_tol;
end
validateattributes(signature.panel_count_candidates,{'numeric'}, ...
    {'vector','integer','positive'},mfilename, ...
    'opt.interpolation_panel_count_candidates');
validateattributes(signature.node_candidates,{'numeric'}, ...
    {'vector','integer','>=',3},mfilename, ...
    'opt.interpolation_node_candidates');
if any(diff(signature.panel_count_candidates)<=0) || ...
        any(diff(signature.node_candidates)<=0)
    error('buildLaplaceCmapInterpolationSignature:UnsortedCandidates', ...
        'Panel and node candidates must be strictly increasing.');
end
if any(mod(signature.node_candidates,2)~=1)
    error('buildLaplaceCmapInterpolationSignature:EvenNodeCount', ...
        'Interpolation node candidates must be odd.');
end
validateattributes(signature.validation_nodes,{'numeric'}, ...
    {'scalar','integer','>',max(signature.node_candidates)},mfilename, ...
    'opt.interpolation_validation_nodes');
end
