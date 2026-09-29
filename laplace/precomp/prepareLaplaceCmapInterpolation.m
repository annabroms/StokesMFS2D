function opt = prepareLaplaceCmapInterpolation(opt,problem)
%PREPARELAPLACECMAPINTERPOLATION Attach a compatible saved model to opt.
%
% Interactive cache misses ask for approval before training. Batch and
% headless runs never train implicitly: call buildLaplaceCmapInterpolator
% explicitly first if a model must be generated noninteractively.

if nargin < 2
    problem = 'capacitance';
end
problem = resolveLaplaceInterpolationProblem(problem,mfilename);
[opt,mode] = configureLaplaceCmapInterpolation(opt,problem);
if strcmp(mode,'none')
    return
end

signature = buildLaplaceCmapInterpolationSignature(opt,true,problem);
signature_id = laplaceInterpolationSignatureId(signature);
if isfield(opt,'interpolation_model') && ...
        ~isempty(opt.interpolation_model)
    model = validatePreparedModel( ...
        opt.interpolation_model,signature,mode,problem);
    opt.interpolation_model = model;
    reportLaplaceCmapInterpolator(model,'provided');
    return
end

model_file = getOptField(opt,'interpolation_model_file','');
if isstring(model_file)
    if ~isscalar(model_file)
        error('prepareLaplaceCmapInterpolation:BadModelFile', ...
            'opt.interpolation_model_file must be scalar text.');
    end
    model_file = char(model_file);
end
use_default_model_file = isempty(model_file);
if use_default_model_file
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    model_dir = fullfile(repo_root,'data','laplace_cmap_interpolation');
    model_file = fullfile(model_dir,sprintf( ...
        'laplace_%s_%s_%s.mat',problem,mode,signature_id));
elseif ~ischar(model_file)
    error('prepareLaplaceCmapInterpolation:BadModelFile', ...
        'opt.interpolation_model_file must be text.');
end
if use_default_model_file && strcmp(problem,'capacitance') && ...
        ~isfile(model_file)
    legacy_signature = legacyCapacitanceSignature(signature);
    legacy_signature_id = laplaceInterpolationSignatureId( ...
        legacy_signature);
    legacy_model_file = fullfile(model_dir,sprintf( ...
        'laplace_capacitance_%s_%s.mat',mode,legacy_signature_id));
    if isfile(legacy_model_file)
        model_file = legacy_model_file;
    end
end
model_file = char(model_file);
opt.interpolation_model_file = model_file;

if isfile(model_file)
    loaded = load(model_file,'model');
    if ~isfield(loaded,'model')
        error('prepareLaplaceCmapInterpolation:BadModelFile', ...
            'The file %s does not contain a struct named model.',model_file);
    end
    model = validatePreparedModel(loaded.model,signature,mode,problem);
    model.model_file = model_file;
    opt.interpolation_model = model;
    markdown_file = replaceExtension(model_file,'.md');
    if ~isfile(markdown_file)
        writeLaplaceCmapInterpolationSummary(model,model_file);
    end
    reportLaplaceCmapInterpolator(model,'loaded');
    return
end

fprintf(2,['\nNo interpolation data exists yet for this model:\n', ...
    '  problem: %s\n  mode: %s\n  C tolerance: %.3e\n', ...
    '  voltage/charge tolerance: %.3e\n', ...
    '  N_cmap/N_f/N_peanut: %d/%d/%d\n', ...
    '  gap range: [%.6g, %.6g]\n  requested file: %s\n'], ...
    problem,mode,opt.interpolation_tol,opt.volt_charge_interp_tol, ...
    opt.N_cmap,opt.N_f,opt.N_peanut,opt.smallest_delta, ...
    opt.delta_pair,model_file);
if ~usejava('desktop')
    error('prepareLaplaceCmapInterpolation:MissingModelBatch', ...
        ['Training is never started implicitly in batch/headless mode. ', ...
         'Run prepareLaplaceCmapInterpolation interactively and approve ', ...
         'training, or explicitly build and save the model first.']);
end

answer = input('Train this interpolation model now? [y/N] ','s');
if ~strcmpi(strtrim(answer),'y') && ~strcmpi(strtrim(answer),'yes')
    error('prepareLaplaceCmapInterpolation:TrainingDeclined', ...
        'Interpolation model training was not approved.');
end

model = buildLaplaceCmapInterpolator(opt,problem);
model.signature = signature;
model.signature_id = signature_id;
model.model_file = model_file;
model_dir = fileparts(model_file);
if ~isempty(model_dir) && ~isfolder(model_dir)
    mkdir(model_dir);
end
save(model_file,'model','-v7.3');
writeLaplaceCmapInterpolationSummary(model,model_file);
opt.interpolation_model = model;
reportLaplaceCmapInterpolator(model,'trained');
end

function model = validatePreparedModel(model,signature,mode,problem)
valid_header = isstruct(model) && isscalar(model) && ...
    isfield(model,'version') && model.version == 3 && ...
    isfield(model,'kind') && isfield(model,'mode') && ...
    strcmp(model.mode,mode) && isfield(model,'signature');
if valid_header
    try
        model_problem = getLaplaceCmapModelProblem(model);
    catch
        model_problem = '';
    end
else
    model_problem = '';
end
if ~valid_header || ~strcmp(model_problem,problem) || ...
        ~isLaplaceCmapSignatureCompatible( ...
        model.signature,signature,problem)
    error('prepareLaplaceCmapInterpolation:IncompatibleModel', ...
        ['The saved/provided interpolation model is incompatible with ', ...
         'the requested problem, mode, tolerances, or construction ', ...
         'parameters.']);
end
if strcmp(problem,'elastance') && ~strcmp(model.mode,'full')
    error('prepareLaplaceCmapInterpolation:ElastanceFullOnly', ...
        'Laplace elastance interpolation models must use mode ''full''.');
end
model.problem = problem;
if ~isfield(model,'volt_charge_action_tolerance') || ...
        isempty(model.volt_charge_action_tolerance)
    model.volt_charge_action_tolerance = ...
        getLaplaceCmapModelVoltChargeTolerance(model);
end
end

function signature = legacyCapacitanceSignature(signature)
tolerance = signature.volt_charge_interp_tol;
signature = rmfield(signature,{'problem','project_charge', ...
    'volt_charge_interp_tol'});
signature.charge_interpolation_tol = tolerance;
end

function output = replaceExtension(filename,new_extension)
[folder,name] = fileparts(filename);
output = fullfile(folder,[name new_extension]);
end
