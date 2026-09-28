function pair_cache = buildLaplaceInterpolatedPairCache(pair_cache,q,pairs, ...
        rbase_in_cmap,rbase_in_f,rout_base_f,rbase_in_solver, ...
        rout_base_c,opt)
%BUILDLAPLACEINTERPOLATEDPAIRCACHE Evaluate a canonical model at pair gaps.

mode = getLaplaceInterpolationMode(opt);
if strcmp(mode,'none')
    error('buildLaplaceInterpolatedPairCache:InterpolationDisabled', ...
        ['This builder requires full, reduced, or reduced_noconst ', ...
         'interpolation.']);
end

timer = tic;
reuse_model = isfield(opt,'interpolation_model') && ...
    ~isempty(opt.interpolation_model);
if reuse_model
    model = opt.interpolation_model;
    validateModel(model,mode,opt,rbase_in_cmap,rbase_in_f,rout_base_f);
else
    error('buildLaplaceInterpolatedPairCache:ModelNotPrepared', ...
        ['Interpolated solves require opt.interpolation_model. Call ', ...
         'opt = prepareLaplaceCmapInterpolation(opt) before the solve.']);
end
model_build_time = toc(timer);

[group_id,group_sep,rep_rows] = groupLaplacePairSeparations( ...
    pair_cache.meta,opt.shared_sep_tol);
n_groups = numel(group_sep);
groups = repmat(emptyInterpolatedGroup(),n_groups,1);
R = opt.rad;
for group_index = 1:n_groups
    delta = group_sep(group_index)-2*R;
    payload = evaluateLaplaceCmapInterpolator(model,delta);
    group = emptyInterpolatedGroup();
    group.group_id = group_index;
    group.sep = group_sep(group_index);
    group.q_pair = [-group.sep/2;group.sep/2];
    group.Cmap = payload.Cmap;
    group.Cmap_QV = payload.Cmap_QV;
    group.rep_pair = pairs(rep_rows(group_index),:);
    group.map_kind = mode;
    group.interp_panel = payload.panel_index;
    group.interp_B = payload.B;
    groups(group_index) = group;
end

pair_cache.enabled = true;
pair_cache.interpolation_mode = mode;
pair_cache.interpolator = model;
pair_cache.groups = groups;
pair_cache.n_groups = n_groups;
pair_cache.group_id = group_id;
pair_cache.group_sep = group_sep;
pair_cache.representative_rows = rep_rows;
pair_cache.stats.branch = ['interpolation_' mode];
pair_cache.stats.n_groups = n_groups;
pair_cache.stats.interpolation_mode = mode;
pair_cache.stats.interpolation_model_reused = reuse_model;
pair_cache.stats.interpolation_model_build_time = model_build_time;
pair_cache.stats.interpolation_exact_snapshots = model.n_exact_snapshots;

for row = 1:size(pairs,1)
    gid = group_id(row);
    pair_cache.meta(row).group_id = gid;
    pair_cache.meta(row).sep = group_sep(gid);
    [Ucross_actual,Ec_actual,Lr_actual] = ...
        buildLaplaceActualPairCollocFactors( ...
        pair_cache.meta(row),q,rbase_in_solver,rout_base_c);
    pair_cache.meta(row).Ucross_colloc_actual = Ucross_actual;
    pair_cache.meta(row).Ec_colloc_actual = Ec_actual;
    pair_cache.meta(row).Lr_colloc_actual = Lr_actual;
end
end

function validateModel(model,mode,opt,rbase_in_c,rbase_in_f,rout_base_f)
if ~isstruct(model) || ~isfield(model,'version') || model.version ~= 3 || ...
        ~isfield(model,'kind') || ...
        ~strcmp(model.kind,'laplace_capacitance_cmap_alpha') || ...
        ~strcmp(model.mode,mode)
    error('buildLaplaceInterpolatedPairCache:ModelModeMismatch', ...
        ['opt.interpolation_model is not a version-3 model for mode ', ...
        '"%s". Rebuild it with buildLaplaceCmapInterpolator.'],mode);
end
signature = buildLaplaceCmapInterpolationSignature(opt,true);
if ~isfield(model,'signature') || ~isequaln(model.signature,signature)
    error('buildLaplaceInterpolatedPairCache:ModelSignatureMismatch', ...
        ['opt.interpolation_model was trained for different tolerances ', ...
         'or canonical construction parameters.']);
end
checks = {model.rbase_in_c,rbase_in_c; ...
    model.rbase_in_f,rbase_in_f; model.rout_base_f,rout_base_f};
for k = 1:size(checks,1)
    expected = checks{k,1}(:);
    actual = checks{k,2}(:);
    if numel(expected) ~= numel(actual) || ...
            norm(expected-actual,inf) > 1e-12*max(1,norm(expected,inf))
        error('buildLaplaceInterpolatedPairCache:ModelGridMismatch', ...
            'opt.interpolation_model was built for different canonical grids.');
    end
end
if abs(model.R-opt.rad) > 1e-12*max(1,opt.rad) || ...
        abs(model.delta_min-opt.smallest_delta) > ...
        1e-12*max(1,opt.smallest_delta) || ...
        abs(model.delta_max-opt.delta_pair) > ...
        1e-12*max(1,opt.delta_pair)
    error('buildLaplaceInterpolatedPairCache:ModelRangeMismatch', ...
        'opt.interpolation_model was built for a different radius or gap range.');
end
end

function group = emptyInterpolatedGroup()
group = struct('group_id',[],'sep',[],'q_pair',[], ...
    'rimage_canon',{{zeros(0,1);zeros(0,1)}}, ...
    'refine_canon',{{zeros(0,1);zeros(0,1)}}, ...
    'Upf',[],'Ypf',[],'DC',[],'YC',[],'Cmap',[], ...
    'Cmap_QV',[],'nout_f',0,'nsrc_f',[0 0],'ntar_f',[0 0], ...
    'rep_pair',[],'map_kind','','interp_panel',[],'interp_B',[]);
end
