%% Tikhonov-regularized Cmap smoothness in the two-body gap
% This is the Tikhonov counterpart of
% sep24_capacitance_cmap_delta_smoothness.m. Two equal disks of radius R
% have centers 0 and 2R+delta. For every gap, the canonical Laplace
% capacitance Cmap is rebuilt with
%
%   opt.use_tikhonov = true,
%   opt.tikhonov_tol = lambda/sigma_max.
%
% The fine-pair and peanut pseudoinverses therefore use the continuous
% filter sigma/(sigma^2+lambda^2), rather than changing retained rank when
% a singular value crosses a TSVD cutoff. Cmap compression is disabled;
% opt.compress_cmap remains a separate truncated-SVD approximation.
%
% By default the ellipse grid is allowed to vary with delta. This retains
% both the main ellipse-segment on/off transition and any smaller source
% count changes caused by individual ellipse nodes crossing the proxy
% filter. The diagnostics distinguish:
%
%   1. the main ellipse on/off transition;
%   2. every ellipse source-count change; and
%   3. any remaining, unexplained nonsmoothness.
%
% A centered interpolation defect in uniformly spaced log10(delta/R),
%   C_k - (C_{k-1}+C_{k+1})/2,
% is used in addition to neighboring steps. It is O(h^2) on a smooth
% trace but becomes large at a jump. Tikhonov effective degrees of freedom
% sum sigma^2/(sigma^2+lambda^2) replace the integer rank plots from sep24.
%
% Anna Broms, Sep 25, 2026

close all;
clear; 

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s ===\n\n',mfilename);

%% Configuration
R = 2;
P = 2;
if ~exist('n_delta','var') || isempty(n_delta)
    n_delta = 300;
end
if ~exist('tikhonov_tol','var') || isempty(tikhonov_tol)
    % Largest value accepted by every solve in the Sep 25 parameter sweep.
    tikhonov_tol = 1e-11;
end
if ~exist('ellipse_constant','var') || isempty(ellipse_constant)
    ellipse_constant = true; % same number of ellipse nodes everywhere, the same placement.
end
if ~exist('plotfig','var') || isempty(plotfig)
    plotfig = true;
end
if ~exist('save_results','var') || isempty(save_results)
    save_results = false;
end

opt = getLaplace2Dparams(P,R);
opt.cmap = true;
opt.compress_cmap = false;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.show_counter = 0;
opt.delta_pair = 0.5*R;
opt.N_f = 150;
opt.use_tikhonov = true;
opt.tikhonov_tol = tikhonov_tol;
opt.ellipse_constant = ellipse_constant;

delta_min = 1e-3*R;
delta_max = 0.3*R;
delta_vals = logspace(log10(delta_min),log10(delta_max),n_delta);
if opt.ellipse_constant
    opt.smallest_delta = delta_min;
end

delta_star = (R-opt.Rp_f)^2/opt.Rp_f;
fprintf(['R=%.3g, N_c=%d, N_f=%d, N_peanut=%d\n', ...
    'delta/R in [%.3e, %.3e], %d log-spaced samples\n', ...
    'use_tikhonov=%d, tikhonov_tol=lambda/sigma_max=%.3e\n', ...
    'ellipse_constant=%d, predicted ellipse switch delta/R=%.5f\n\n'], ...
    R,opt.N_c,opt.N_f,opt.N_peanut,delta_min/R,delta_max/R,n_delta, ...
    opt.use_tikhonov,opt.tikhonov_tol,opt.ellipse_constant,delta_star/R);

%% Fixed coarse/fine circular grids
N_c = opt.N_c;
N_f = opt.N_f;
nout_c = ceil(opt.a_c*N_c);
nout_f = ceil(opt.a_f*N_f);

tout_c = (0:nout_c-1)'*(2*pi/nout_c);
rbase_out_c = R*exp(1i*tout_c);
tin_c = (0:N_c-1)'*(2*pi/N_c);
rbase_in_c = opt.Rp_c*exp(1i*tin_c);
tin_f = (0:N_f-1)'*(2*pi/N_f);
rbase_in_f = opt.Rp_f*exp(1i*tin_f);
tout_f = (0:nout_f-1)'*(2*pi/nout_f);
rout_base_f = R*exp(1i*tout_f);

%% Storage and representative entries
N = 2*N_c;
idx_facing = round(N_c/2)+1;
track_idx = [ ...
    1,1; ...
    1,N_c+idx_facing; ...
    N_c+idx_facing,N_c+idx_facing; ...
    idx_facing,N_c+1; ...
    1,N_c+1; ...
    idx_facing,N_c+idx_facing];
n_track = size(track_idx,1);
track_labels = arrayfun(@(k) sprintf('C(%d,%d)', ...
    track_idx(k,1),track_idx(k,2)),1:n_track,'UniformOutput',false);

C_all = nan(N,N,n_delta);
entries = nan(n_delta,n_track);
fro_norm = nan(n_delta,1);
max_abs = nan(n_delta,1);
n_image = zeros(n_delta,2);
effective_dof_fine = nan(n_delta,1);
effective_dof_peanut = nan(n_delta,1);
n_modes_fine = nan(n_delta,1);
n_modes_peanut = nan(n_delta,1);

fprintf('Building Tikhonov-regularized Cmap at every gap...\n');
for k = 1:n_delta
    delta = delta_vals(k);
    q = [0; 2*R+delta];

    [~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
    assert(isequal(pairs,[1 2]),'Expected the single close pair (1,2).');
    n_image(k,:) = [numel(rimage_vec{1,2}),numel(rimage_vec{2,1})];

    [~,~,~,~,Cmap] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f, ...
        rout_base_f,rbase_out_c,rimage_vec,refine,pairs,opt,rbase_in_c);
    C = Cmap{1,2};
    C_all(:,:,k) = C;
    fro_norm(k) = norm(C,'fro');
    max_abs(k) = max(abs(C(:)));
    for m = 1:n_track
        entries(k,m) = C(track_idx(m,1),track_idx(m,2));
    end

    % Replicate the two matrices regularized in the production Cmap build
    % solely to report their continuous effective degrees of freedom.
    rin_pair_f = [q(1)+rbase_in_f; rimage_vec{1,2}; ...
        q(2)+rbase_in_f; rimage_vec{2,1}];
    rout_f = [q(1)+rout_base_f; refine{1,2}; ...
        q(2)+rout_base_f; refine{2,1}];
    s_fine = svd(lapSLPmat(rin_pair_f,rout_f));
    effective_dof_fine(k) = tikhonov_effective_dof(s_fine,tikhonov_tol);
    n_modes_fine(k) = numel(s_fine);

    rin_pair_c = [q(1)+rbase_in_c; q(2)+rbase_in_c];
    rout_peanut = createPeanut(q(1),q(2),opt.N_peanut,0,R);
    s_peanut = svd(lapSLPmat(rin_pair_c,rout_peanut));
    effective_dof_peanut(k) = tikhonov_effective_dof( ...
        s_peanut,tikhonov_tol);
    n_modes_peanut(k) = numel(s_peanut);
end
fprintf('Done.\n\n');

%% Locate ellipse discretisation changes
ellipse_switch_idx = find(any(diff(n_image,1,1)~=0,2));
is_enhanced = any(n_image>0,2);
activation_switch_idx = find(diff(is_enhanced)~=0);
fprintf('Ellipse source-count changes: %d\n',numel(ellipse_switch_idx));
fprintf('Whole ellipse on/off transitions: %d\n',numel(activation_switch_idx));
if ~isempty(ellipse_switch_idx)
    fprintf('All ellipse switch locations delta/R: %s\n', ...
        mat2str(round(delta_vals(ellipse_switch_idx+1)/R,5)));
end
if ~isempty(activation_switch_idx)
    fprintf('Main on/off location delta/R: %s\n', ...
        mat2str(round(delta_vals(activation_switch_idx+1)/R,5)));
end
fprintf('\n');

%% Neighbor steps and centered interpolation defects for all Cmap entries
V = reshape(C_all,N*N,n_delta);
entry_scale = max(abs(V),[],2);
entry_scale(entry_scale==0) = 1;

step_fro = nan(n_delta-1,1);
step_entry_max = nan(n_delta-1,1);
for k = 1:n_delta-1
    dC = C_all(:,:,k+1)-C_all(:,:,k);
    scale = max([norm(C_all(:,:,k),'fro'),norm(C_all(:,:,k+1),'fro'),eps]);
    step_fro(k) = norm(dC,'fro')/scale;
    step_entry_max(k) = max(abs(V(:,k+1)-V(:,k))./entry_scale);
end

defect_fro = nan(n_delta,1);
defect_entry_max = nan(n_delta,1);
for k = 2:n_delta-1
    defect = C_all(:,:,k)-0.5*(C_all(:,:,k-1)+C_all(:,:,k+1));
    defect_fro(k) = norm(defect,'fro')/max(norm(C_all(:,:,k),'fro'),eps);
    entry_defect = V(:,k)-0.5*(V(:,k-1)+V(:,k+1));
    defect_entry_max(k) = max(abs(entry_defect)./entry_scale);
end

center_idx = (2:n_delta-1)';
near_activation = near_step_switch(center_idx,activation_switch_idx,1);
near_any_ellipse = near_step_switch(center_idx,ellipse_switch_idx,1);
valid = isfinite(defect_entry_max(center_idx));

[largest_all,idx_all] = max_with_index(defect_entry_max,valid,center_idx);
[largest_off_activation,idx_off_activation] = max_with_index( ...
    defect_entry_max,valid & ~near_activation,center_idx);
[largest_unexplained,idx_unexplained] = max_with_index( ...
    defect_entry_max,valid & ~near_any_ellipse,center_idx);

background_values = defect_entry_max(center_idx(valid & ~near_any_ellipse));
background_median = median(background_values,'omitnan');
background_p99 = percentile(background_values,99);
unexplained_threshold = max(50*background_median,5*background_p99);
unexplained_centers = center_idx(valid & ~near_any_ellipse & ...
    defect_entry_max(center_idx)>unexplained_threshold);

fprintf('Largest max-entry interpolation defect over all gaps: %.3e at delta/R %.6f\n', ...
    largest_all,delta_vals(idx_all)/R);
fprintf(['Largest defect away from only the main ellipse on/off transition: ', ...
    '%.3e at delta/R %.6f\n'], ...
    largest_off_activation,delta_vals(idx_off_activation)/R);
fprintf(['Largest defect away from every ellipse source-count change: ', ...
    '%.3e at delta/R %.6f\n'], ...
    largest_unexplained,delta_vals(idx_unexplained)/R);
fprintf('Away-switch background: median %.3e, 99th percentile %.3e\n', ...
    background_median,background_p99);
fprintf('Unexplained outliers above %.3e: %d\n\n', ...
    unexplained_threshold,numel(unexplained_centers));

% Print the strongest local defects and whether an ellipse change explains
% them. This makes the conclusion reproducible without inspecting figures.
[~,order] = sort(defect_entry_max(center_idx),'descend','MissingPlacement','last');
n_report = min(12,numel(order));
fprintf('%12s %14s %12s %12s\n', ...
    'delta/R','entry defect','main switch','any ellipse');
for j = 1:n_report
    k = center_idx(order(j));
    fprintf('%12.6f %14.3e %12s %12s\n',delta_vals(k)/R, ...
        defect_entry_max(k),mat2str(near_step_switch(k,activation_switch_idx,1)), ...
        mat2str(near_step_switch(k,ellipse_switch_idx,1)));
end
fprintf('\n');

%% Figures
if plotfig
    figure('Name','Tikhonov Cmap entries vs gap','Color','w');
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    nexttile;
    semilogx(delta_vals/R,entries,'LineWidth',1.0);
    hold on;
    add_switch_lines(delta_vals,R,ellipse_switch_idx,activation_switch_idx);
    xlabel('\delta/R'); ylabel('Cmap entry');
    title(sprintf('Tikhonov Cmap entries, \eta=%.1e',tikhonov_tol));
    legend(track_labels,'Location','best'); grid on;
    nexttile;
    semilogx(delta_vals/R,fro_norm,'DisplayName','||Cmap||_F');
    hold on;
    semilogx(delta_vals/R,max_abs,'DisplayName','max |Cmap|');
    add_switch_lines(delta_vals,R,ellipse_switch_idx,activation_switch_idx);
    xlabel('\delta/R'); ylabel('aggregate magnitude');
    legend('Location','best'); grid on;

    figure('Name','Ellipse counts and Tikhonov effective DOF','Color','w');
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    nexttile;
    stairs(delta_vals/R,n_image(:,1),'DisplayName','body 1'); hold on;
    stairs(delta_vals/R,n_image(:,2),'DisplayName','body 2');
    set(gca,'XScale','log'); xlabel('\delta/R');
    ylabel('ellipse source count'); legend('Location','best'); grid on;
    nexttile;
    semilogx(delta_vals/R,effective_dof_fine,'DisplayName','fine pair'); hold on;
    semilogx(delta_vals/R,effective_dof_peanut,'DisplayName','peanut');
    add_switch_lines(delta_vals,R,ellipse_switch_idx,activation_switch_idx);
    xlabel('\delta/R'); ylabel('effective degrees of freedom');
    legend('Location','best'); grid on;

    figure('Name','Tikhonov Cmap jump diagnostics','Color','w');
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    nexttile;
    loglog(delta_vals(2:end)/R,step_entry_max,'DisplayName','max entry step');
    hold on;
    loglog(delta_vals(2:end)/R,step_fro,'DisplayName','relative Frobenius step');
    add_switch_lines(delta_vals,R,ellipse_switch_idx,activation_switch_idx);
    xlabel('\delta/R'); ylabel('neighboring step');
    legend('Location','best'); grid on;
    nexttile;
    loglog(delta_vals/R,defect_entry_max,'DisplayName','max entry defect');
    hold on;
    loglog(delta_vals/R,defect_fro,'DisplayName','relative Frobenius defect');
    add_switch_lines(delta_vals,R,ellipse_switch_idx,activation_switch_idx);
    xlabel('\delta/R'); ylabel('centered interpolation defect');
    legend('Location','best'); grid on;

    V_range = max(V,[],2)-min(V,[],2);
    V_range(V_range==0) = 1;
    V_norm = (V-min(V,[],2))./V_range;
    figure('Name','All Tikhonov Cmap entries vs gap','Color','w');
    imagesc(log10(delta_vals/R),1:size(V,1),V_norm);
    set(gca,'YDir','normal'); colorbar;
    xlabel('log_{10}(\delta/R)');
    ylabel('flattened Cmap entry index');
    title('Every Tikhonov Cmap entry, individually rescaled to [0,1]');
end

result = struct();
result.options = opt;
result.delta_vals = delta_vals;
result.tikhonov_tol = tikhonov_tol;
result.ellipse_constant = opt.ellipse_constant;
result.delta_star = delta_star;
result.n_image = n_image;
result.ellipse_switch_idx = ellipse_switch_idx;
result.activation_switch_idx = activation_switch_idx;
result.effective_dof_fine = effective_dof_fine;
result.effective_dof_peanut = effective_dof_peanut;
result.n_modes_fine = n_modes_fine;
result.n_modes_peanut = n_modes_peanut;
result.entries = entries;
result.track_idx = track_idx;
result.fro_norm = fro_norm;
result.max_abs = max_abs;
result.step_fro = step_fro;
result.step_entry_max = step_entry_max;
result.defect_fro = defect_fro;
result.defect_entry_max = defect_entry_max;
result.background_median = background_median;
result.background_p99 = background_p99;
result.unexplained_threshold = unexplained_threshold;
result.unexplained_centers = unexplained_centers;
result.largest_all = largest_all;
result.largest_off_activation = largest_off_activation;
result.largest_unexplained = largest_unexplained;

if save_results
    result_file = fullfile(fileparts(mfilename('fullpath')), ...
        'sep25_capacitance_cmap_delta_smoothness_tikhonov_results.mat');
    save(result_file,'result','-v7.3');
    fprintf('Saved %s\n',result_file);
end

fprintf('Results are available in workspace variable result.\n');

function dof = tikhonov_effective_dof(s,eta)
if isempty(s) || max(s)==0
    dof = 0;
    return
end
s = s/max(s);
dof = sum(s.^2./(s.^2+eta^2));
end

function mask = near_step_switch(center_idx,switch_idx,radius)
center_idx = center_idx(:);
switch_idx = switch_idx(:)';
if isempty(switch_idx)
    mask = false(size(center_idx));
    return
end
% A switch index s denotes the step s -> s+1. A centered defect at either
% adjacent sample can contain that step; radius expands this neighborhood.
distance = min(abs(center_idx-switch_idx),abs(center_idx-(switch_idx+1)));
mask = any(distance<=radius,2);
end

function [value,index] = max_with_index(values,mask,candidate_idx)
idx = candidate_idx(mask);
if isempty(idx)
    value = NaN;
    index = NaN;
    return
end
[value,local] = max(values(idx));
index = idx(local);
end

function value = percentile(values,p)
values = sort(values(isfinite(values)));
if isempty(values)
    value = NaN;
    return
end
position = 1+(numel(values)-1)*p/100;
lo = floor(position);
hi = ceil(position);
if lo==hi
    value = values(lo);
else
    value = values(lo)+(position-lo)*(values(hi)-values(lo));
end
end

function add_switch_lines(delta_vals,R,ellipse_switch_idx,activation_switch_idx)
for s = ellipse_switch_idx(:)'
    xline(delta_vals(s+1)/R,':','Color',[0.65 0.65 0.65], ...
        'HandleVisibility','off');
end
for s = activation_switch_idx(:)'
    xline(delta_vals(s+1)/R,'--','Color',[0.85 0.20 0.10], ...
        'LineWidth',1.2,'HandleVisibility','off');
end
end
