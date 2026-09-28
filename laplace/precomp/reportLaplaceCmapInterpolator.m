function reportLaplaceCmapInterpolator(model,source)
%REPORTLAPLACECMAPINTERPOLATOR Print the trained/loaded model summary.

if nargin < 2
    source = 'model';
end
fprintf('\n=== Laplace capacitance interpolation model (%s) ===\n',source);
fprintf(['mode=%s, C tolerance %.3e, C_Q tolerance %.3e, ', ...
    'panels=%d, exact snapshots=%d\n'],model.mode, ...
    model.action_tolerance,model.charge_action_tolerance, ...
    model.n_panels,model.n_exact_snapshots);
fprintf('C nodes by panel:   %s\n',mat2str(model.q_by_panel));
fprintf('C_Q nodes by panel: %s\n',mat2str(model.q_charge_by_panel));
if any(strcmp(model.mode,{'reduced','reduced_noconst'}))
    fprintf('spatial ranks:      %s\n',mat2str(model.spatial_ranks));
else
    fprintf('spatial ranks:      not used by the full model\n');
end
fprintf('panel counts by base range: %s\n', ...
    mat2str(model.panel_counts_by_base));
fprintf('certified max C error %.3e; max C_Q error %.3e\n', ...
    model.training.max_C_error,model.training.max_CQ_error);
if isfield(model,'model_file') && ~isempty(model.model_file)
    fprintf('MAT file: %s\n',model.model_file);
    [folder,name] = fileparts(model.model_file);
    fprintf('summary:  %s\n',fullfile(folder,[name '.md']));
end
fprintf('==========================================================\n\n');
end
