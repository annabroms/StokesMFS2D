function markdown_file = writeLaplaceCmapInterpolationSummary(model,model_file)
%WRITELAPLACECMAPINTERPOLATIONSUMMARY Write a human-readable model sidecar.

[folder,name] = fileparts(model_file);
markdown_file = fullfile(folder,[name '.md']);
fid = fopen(markdown_file,'w');
if fid < 0
    error('writeLaplaceCmapInterpolationSummary:OpenFailed', ...
        'Could not open %s for writing.',markdown_file);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid,'# Laplace capacitance interpolation model\n\n');
fprintf(fid,'- Mode: `%s`\n',model.mode);
fprintf(fid,'- Signature: `%s`\n',model.signature_id);
fprintf(fid,'- Matrix tolerance: `%.16g`\n',model.action_tolerance);
fprintf(fid,'- Charge-map tolerance: `%.16g`\n', ...
    model.charge_action_tolerance);
fprintf(fid,'- Coordinate: `alpha`\n');
fprintf(fid,'- Gap interval: `[%.16g, %.16g]`\n', ...
    model.delta_min,model.delta_max);
fprintf(fid,'- Panels: `%d` (`%s` by base range)\n', ...
    model.n_panels,strtrim(sprintf('%d ',model.panel_counts_by_base)));
fprintf(fid,'- Exact snapshots used: `%d`\n',model.n_exact_snapshots);
fprintf(fid,'- Maximum certified C error: `%.6e`\n', ...
    model.training.max_C_error);
fprintf(fid,'- Maximum certified C_Q error: `%.6e`\n\n', ...
    model.training.max_CQ_error);

fprintf(fid,'## Discretization and regularization\n\n');
fprintf(fid,'| Parameter | Value |\n|---|---:|\n');
fields = {'rad','N_cmap','N_f','N_peanut','a_f','Rp_c','Rp_f', ...
    'delta_min','delta_max','beta','Nclust','tikhonov_tol'};
for k = 1:numel(fields)
    field = fields{k};
    fprintf(fid,'| `%s` | `%.16g` |\n',field,model.signature.(field));
end

fprintf(fid,'\n## Selected panels\n\n');
fprintf(fid,['| Panel | Base | delta min | delta max | C nodes | ', ...
    'C_Q nodes | Rank | C error | C_Q error |\n']);
fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for k = 1:model.n_panels
    panel = model.panels(k);
    if isempty(panel.spatial_rank)
        rank_value = '--';
    else
        rank_value = sprintf('%d',panel.spatial_rank);
    end
    fprintf(fid,['| %d | %d | %.6e | %.6e | %d | %d | %s | ', ...
        '%.6e | %.6e |\n'],k,panel.base_index,panel.delta_lo, ...
        panel.delta_hi,panel.q,panel.q_charge,rank_value, ...
        panel.certified_C_error,panel.certified_CQ_error);
end

fprintf(fid,'\n## MATLAB usage\n\n```matlab\n');
fprintf(fid,'opt.use_interpolation = ''%s'';\n',model.mode);
fprintf(fid,'loaded = load(''%s'',''model'');\n',strrep(model_file,'''',''''''));
fprintf(fid,'opt.interpolation_model = loaded.model;\n```\n');
clear cleanup
end
