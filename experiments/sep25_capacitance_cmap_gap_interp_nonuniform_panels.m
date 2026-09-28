%% SEP25_CAPACITANCE_CMAP_GAP_INTERP_NONUNIFORM_PANELS
% Count-aware nonuniform interpolation panels for the ellipse-enhanced
% Tikhonov-regularised Laplace capacitance pair map.
%
% Uniform panels can straddle a gap at which one pair of ellipse sources
% crosses the fine proxy circle.  The matrix construction then changes its
% discrete source set inside one polynomial panel.  This experiment first
% locates every ellipse-source-count transition and uses those locations as
% nonuniform numerical panel boundaries.  Wider constant-count intervals
% are then subdivided until their log-gap width is at most 0.08.  The script
% chooses the smallest q=3,5,9,17 Chebyshev--Lobatto interpolant that meets
% a 1e-6 relative spectral/action-error target on each resulting panel.
%
% A relative margin is removed on either side of each count transition so
% that no Lobatto endpoint sits on the branch condition.  The omitted
% brackets are only transition_margin_s wide in log-gap and can be handled
% by an exact Cmap build (or a separately defined one-sided convention) in
% an online implementation.  They are not physical regimes.
%
% Figures are displayed directly and also saved under data/.
%
% Anna Broms, Sep 25, 2026

close all;
clearvars;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;
addpath(fullfile(repo_root,'experiments','cap_interp'));

fprintf('=== sep25_capacitance_cmap_gap_interp_nonuniform_panels ===\n\n');

%% Configuration
R = 2;
opt = getLaplace2Dparams(2,R);
opt.cmap = 1;
opt.compress_cmap = false;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.show_counter = 0;
opt.delta_pair = 0.2*R;
opt.ellipse_constant = false;
opt.use_tikhonov = true;
opt.tikhonov_tol = 1e-11;
opt.Nclust = 150;

action_tolerance = 1e-6;
q_candidates = [3 5 9];
q_validation = 33;
n_random_per_panel = 4;
rng_seed = 20260925;
n_transition_scan = 4001;
transition_margin_s = 1e-9;
physical_transition_margin = 1e-8;
max_panel_width_s = 0.08;

delta_min = 1e-3*R;
delta_star = (R-opt.Rp_f)^2/opt.Rp_f;
delta_hi = delta_star*(1-physical_transition_margin);
to_param = @(delta) log(delta/R);
to_delta = @(s) R*exp(s);
s_lo = to_param(delta_min);
s_hi = to_param(delta_hi);

cache_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_tikhonov_cache.mat');
results_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_nonuniform_panels_results.mat');
figure_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_nonuniform_panels.png');

%% Fixed discretisation grids
N_c = opt.N_c;
N_f = opt.N_f;
nout_c = ceil(opt.a_c*N_c);
tout_c = linspace(0,2*pi,nout_c+1)';
tout_c = tout_c(1:end-1);
grids.rbase_out_c = R*(cos(tout_c)+1i*sin(tout_c));

tin_c = linspace(0,2*pi,N_c+1)';
tin_c = tin_c(1:end-1);
grids.rbase_in_c = opt.Rp_c*(cos(tin_c)+1i*sin(tin_c));

tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
grids.rbase_in_f = opt.Rp_f*(cos(tin_f)+1i*sin(tin_f));

nout_f = ceil(opt.a_f*N_f);
tout_f = linspace(0,2*pi,nout_f+1)';
tout_f = tout_f(1:end-1);
grids.rout_base_f = R*(cos(tout_f)+1i*sin(tout_f));

%% Locate every discrete ellipse-source-count transition
[transition_s,scan_s,scan_count] = findCountTransitions( ...
    s_lo,s_hi,n_transition_scan,R,opt,to_delta);
count_panel_lo = [s_lo; transition_s+transition_margin_s];
count_panel_hi = [transition_s-transition_margin_s; s_hi];
assert(all(count_panel_hi>count_panel_lo), ...
    'Transition margin produced an empty count-aware interval.');

panel_lo = zeros(0,1);
panel_hi = zeros(0,1);
for k = 1:numel(count_panel_lo)
    n_piece = ceil((count_panel_hi(k)-count_panel_lo(k))/max_panel_width_s);
    edges = linspace(count_panel_lo(k),count_panel_hi(k),n_piece+1)';
    panel_lo = [panel_lo; edges(1:end-1)]; %#ok<AGROW>
    panel_hi = [panel_hi; edges(2:end)]; %#ok<AGROW>
end
n_panels = numel(panel_lo);

fprintf(['Nclust=%d gives %d source-count transitions and %d nonuniform ', ...
    'panels after imposing max width %.3g in log(delta/R).\n'], ...
    opt.Nclust,numel(transition_s),n_panels,max_panel_width_s);
fprintf('transition delta/R = %s\n\n',mat2str(exp(transition_s),7));

%% Plan all exact validation and random snapshots
rng(rng_seed);
plans = repmat(struct('s_lo',[],'s_hi',[],'validation_s',[], ...
    'random_s',[]),n_panels,1);
all_delta = zeros(0,1);
for ip = 1:n_panels
    plans(ip).s_lo = panel_lo(ip);
    plans(ip).s_hi = panel_hi(ip);
    plans(ip).validation_s = chebLobattoNodes( ...
        panel_lo(ip),panel_hi(ip),q_validation);
    plans(ip).random_s = panel_lo(ip)+(panel_hi(ip)-panel_lo(ip))* ...
        rand(n_random_per_panel,1);
    all_delta = [all_delta; to_delta(plans(ip).validation_s); ...
        to_delta(plans(ip).random_s)]; %#ok<AGROW>
end
delta_pool = mergeNearlyEqual(all_delta,1e-12);

fprintf('Building/reusing %d exact snapshots...\n',numel(delta_pool));
[C_pool,info_pool] = getOrBuildCmapSnapshots( ...
    delta_pool,R,opt,grids,cache_file);

%% Select the smallest passing local q on every nonuniform panel
q_local = nan(1,n_panels);
held_error = nan(1,n_panels);
random_error = nan(1,n_panels);
error_by_q = nan(n_panels,numel(q_candidates));
panel_count = nan(1,n_panels);

for ip = 1:n_panels
    plan = plans(ip);
    C_validation = lookupMatrices( ...
        to_delta(plan.validation_s),delta_pool,C_pool);

    count_at_nodes = ellipseCounts(to_delta(plan.validation_s),R,opt);
    assert(all(count_at_nodes==count_at_nodes(1)), ...
        'A count-aware panel still straddles an ellipse-count transition.');
    panel_count(ip) = count_at_nodes(1);

    for iq = 1:numel(q_candidates)
        q = q_candidates(iq);
        s_train = chebLobattoNodes(plan.s_lo,plan.s_hi,q);
        C_train = lookupMatrices(to_delta(s_train),delta_pool,C_pool);
        held = heldOutMask(plan.validation_s,s_train, ...
            plan.s_hi-plan.s_lo);
        C_hat = evalMatrixChebBary(plan.validation_s(held),s_train, ...
            chebBarycentricWeights(q),C_train);
        error_by_q(ip,iq) = max(relativeActionError( ...
            C_hat,C_validation(:,:,held)));
    end

    iq = find(error_by_q(ip,:) <= action_tolerance,1,'first');
    if isempty(iq)
        continue
    end
    q_local(ip) = q_candidates(iq);
    held_error(ip) = error_by_q(ip,iq);

    s_train = chebLobattoNodes(plan.s_lo,plan.s_hi,q_local(ip));
    C_train = lookupMatrices(to_delta(s_train),delta_pool,C_pool);
    C_random = lookupMatrices(to_delta(plan.random_s),delta_pool,C_pool);
    C_hat_random = evalMatrixChebBary(plan.random_s,s_train, ...
        chebBarycentricWeights(q_local(ip)),C_train);
    random_error(ip) = max(relativeActionError(C_hat_random,C_random));
end

pass = all(isfinite(q_local)) && max(random_error) <= action_tolerance;
fprintf('\n%-5s %-21s %6s %7s %11s %11s\n', ...
    'panel','delta/R interval','nodes','q','held-out','random');
for ip = 1:n_panels
    fprintf('%-5d [%.4e, %.4e] %6d %7s %11s %11s\n',ip, ...
        exp(panel_lo(ip)),exp(panel_hi(ip)),panel_count(ip), ...
        numberOrDash(q_local(ip),'%d'), ...
        numberOrDash(held_error(ip),'%.3e'), ...
        numberOrDash(random_error(ip),'%.3e'));
end
fprintf(['\nOverall: %s; q_max=%s; held-out max=%s; random max=%s ', ...
    '(target %.1e).\n'],passText(pass), ...
    numberOrDash(max(q_local,[],'omitnan'),'%d'), ...
    numberOrDash(max(held_error,[],'omitnan'),'%.3e'), ...
    numberOrDash(max(random_error,[],'omitnan'),'%.3e'),action_tolerance);

%% Visible diagnostic figures
f = figure('Name','Nonuniform count-aware Cmap interpolation panels', ...
    'Color','w','Visible','on');
subplot(2,1,1);
semilogx(exp(scan_s),scan_count,'LineWidth',1.1); hold on;
for k = 1:numel(transition_s)
    xline(exp(transition_s(k)),'--');
end
xlabel('\delta/R'); ylabel('retained ellipse sources per body');
title('Count-aware nonuniform panel boundaries'); grid on;

subplot(2,1,2);
centres = exp((panel_lo+panel_hi)/2);
semilogy(centres,held_error,'-o','LineWidth',1.2, ...
    'DisplayName','held-out'); hold on;
semilogy(centres,random_error,'-s','LineWidth',1.2, ...
    'DisplayName','random');
yline(action_tolerance,'--','DisplayName','target');
xlabel('panel-centre \delta/R'); ylabel('relative action error');
title(sprintf('Selected local orders q=%s',mat2str(q_local)));
legend('Location','best'); grid on;

saveas(f,figure_file);
drawnow;

retained_pool = arrayfun(@(x) x.n_image(1),info_pool);
save(results_file,'pass','q_local','held_error','random_error', ...
    'error_by_q','panel_lo','panel_hi','panel_count','transition_s', ...
    'scan_s','scan_count','retained_pool','action_tolerance', ...
    'q_candidates','q_validation','n_random_per_panel','rng_seed', ...
    'transition_margin_s','physical_transition_margin', ...
    'max_panel_width_s','opt','R','-v7.3');
fprintf('Results saved to %s\n',relativePath(results_file,repo_root));
fprintf('Figure saved to %s\n',relativePath(figure_file,repo_root));


%% Local helpers
function [switch_s,scan_s,scan_count] = findCountTransitions( ...
        s_lo,s_hi,n_scan,R,opt,to_delta)
scan_s = linspace(s_lo,s_hi,n_scan)';
scan_count = ellipseCounts(to_delta(scan_s),R,opt);
indices = find(diff(scan_count)~=0);
switch_s = zeros(numel(indices),1);
for k = 1:numel(indices)
    left = scan_s(indices(k));
    right = scan_s(indices(k)+1);
    count_left = scan_count(indices(k));
    for iteration = 1:50
        midpoint = (left+right)/2;
        count_mid = ellipseCounts(to_delta(midpoint),R,opt);
        if count_mid==count_left
            left = midpoint;
        else
            right = midpoint;
        end
    end
    switch_s(k) = (left+right)/2;
end
end

function count = ellipseCounts(delta,R,opt)
delta = delta(:);
count = zeros(size(delta));
for k = 1:numel(delta)
    q_centres = [0;2*R+delta(k)];
    [~,~,~,rimage_vec] = getEnhancedGrid(q_centres,opt);
    count(k) = numel(rimage_vec{1,2});
end
end

function mask = heldOutMask(validation_nodes,training_nodes,scale)
distance = min(abs(validation_nodes(:)-training_nodes(:)'),[],2);
mask = distance > 1e-10*max(1,abs(scale));
end

function error = relativeActionError(C_hat,C_exact)
n_test = size(C_exact,3);
error = zeros(n_test,1);
for k = 1:n_test
    error(k) = norm(C_hat(:,:,k)-C_exact(:,:,k),2)/ ...
        max(norm(C_exact(:,:,k),2),eps);
end
end

function values = mergeNearlyEqual(values,relative_tolerance)
values = sort(values(:));
if isempty(values)
    return
end
keep = true(size(values));
last = values(1);
for k = 2:numel(values)
    if abs(values(k)-last) <= relative_tolerance* ...
            max([1,abs(values(k)),abs(last)])
        keep(k) = false;
    else
        last = values(k);
    end
end
values = values(keep);
end

function C = lookupMatrices(delta_query,delta_pool,C_pool)
delta_query = delta_query(:);
n = size(C_pool,1);
C = zeros(n,n,numel(delta_query));
for k = 1:numel(delta_query)
    [distance,index] = min(abs(delta_pool-delta_query(k)));
    assert(distance <= 1e-10*max(1,abs(delta_query(k))), ...
        'Requested delta is absent from the exact-snapshot pool.');
    C(:,:,k) = C_pool(:,:,index);
end
end

function value = numberOrDash(x,format)
if isfinite(x)
    value = sprintf(format,x);
else
    value = '--';
end
end

function value = passText(pass)
if pass
    value = 'PASS';
else
    value = 'FAIL';
end
end

function p = relativePath(filename,root)
prefix = [root filesep];
if strncmp(filename,prefix,numel(prefix))
    p = filename(numel(prefix)+1:end);
else
    p = filename;
end
end
