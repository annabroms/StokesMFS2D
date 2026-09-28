function opt = prepareLaplaceCmapInterpolation(opt)
%PREPARELAPLACECMAPINTERPOLATION Attach a compatible saved model to opt.
%
% Interactive cache misses ask for approval before training.  Batch and
% headless runs never train implicitly: call buildLaplaceCmapInterpolator
% explicitly first if a model must be generated noninteractively.

[opt,mode] = configureLaplaceCapacitanceInterpolation(opt);
if strcmp(mode,'none')
    return
end

signature = buildLaplaceCmapInterpolationSignature(opt,true);
signature_id = laplaceInterpolationSignatureId(signature);
if isfield(opt,'interpolation_model') && ...
        ~isempty(opt.interpolation_model)
    validatePreparedModel(opt.interpolation_model,signature,mode);
    reportLaplaceCmapInterpolator(opt.interpolation_model,'provided');
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
if isempty(model_file)
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    model_dir = fullfile(repo_root,'data','laplace_cmap_interpolation');
    model_file = fullfile(model_dir,sprintf( ...
        'laplace_capacitance_%s_%s.mat',mode,signature_id));
elseif ~ischar(model_file)
    error('prepareLaplaceCmapInterpolation:BadModelFile', ...
        'opt.interpolation_model_file must be text.');
end
model_file = char(model_file);
opt.interpolation_model_file = model_file;

if isfile(model_file)
    loaded = load(model_file,'model');
    if ~isfield(loaded,'model')
        error('prepareLaplaceCmapInterpolation:BadModelFile', ...
            'The file %s does not contain a struct named model.',model_file);
    end
    validatePreparedModel(loaded.model,signature,mode);
    model = loaded.model;
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
    '  mode: %s\n  C tolerance: %.3e\n', ...
    '  C_Q tolerance: %.3e\n  N_cmap/N_f/N_peanut: %d/%d/%d\n', ...
    '  gap range: [%.6g, %.6g]\n  requested file: %s\n'], ...
    mode,opt.interpolation_tol,opt.charge_interpolation_tol, ...
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

model = buildLaplaceCmapInterpolator(opt);
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

function validatePreparedModel(model,signature,mode)
if ~isstruct(model) || ~isfield(model,'version') || model.version ~= 3 || ...
        ~isfield(model,'kind') || ...
        ~strcmp(model.kind,'laplace_capacitance_cmap_alpha') || ...
        ~isfield(model,'mode') || ~strcmp(model.mode,mode) || ...
        ~isfield(model,'signature') || ...
        ~isequaln(model.signature,signature)
    error('prepareLaplaceCmapInterpolation:IncompatibleModel', ...
        ['The saved/provided interpolation model is incompatible with ', ...
         'the requested mode, tolerances, or construction parameters.']);
end
end

function output = replaceExtension(filename,new_extension)
[folder,name] = fileparts(filename);
output = fullfile(folder,[name new_extension]);
end
