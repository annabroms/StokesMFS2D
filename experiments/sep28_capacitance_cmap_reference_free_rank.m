%% SEP28_CAPACITANCE_CMAP_REFERENCE_FREE_RANK
% Test the panel-local ranks required by a reference-free approximation
%
%   C(alpha) approximately U_p*B_p(alpha)*V_p',
%
% instead of the current centered approximation
%
%   C(alpha) approximately C_ref,p + U_p*B_p(alpha)*V_p'.
%
% The test uses the twelve production alpha panels (four below and eight
% above ellipse activation), Tikhonov regularisation, a constant ellipse
% discretisation with 150 nodes, and no Cmap TSVD compression.  Common
% left/right bases are trained directly on uncentered C snapshots.  Ranks
% are audited at held-out Chebyshev--Lobatto nodes and deterministic
% off-grid alpha values using relative spectral/action error.
%
% Two ranks are reported separately:
%   1. projection-only rank, using exact C at every audit point;
%   2. end-to-end rank after B(alpha) is barycentrically interpolated.
%
% The production exact-snapshot cache is reused when compatible.  Missing
% snapshots fall back to a small experiment-local cache.  Figures are
% displayed directly and also saved under data/.
%
% Anna Broms, Sep 28, 2026

close all;
clearvars;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;
addpath(fullfile(repo_root,'experiments','cap_interp'));

fprintf('=== sep28_capacitance_cmap_reference_free_rank ===\n\n');

%% Configuration
R = 2;
P = 2;
action_tolerance = 1e-6;
spatial_action_tolerance = 2.5e-7;
q_basis = 17;
q_validation = 33;
q_candidates = [5 7 9 17];
audit_fraction = [0.137 0.367 0.643 0.883];
rank_search_stride = 4;
rng_seed = 20260928;

opt = getLaplace2Dparams(P,R);
opt.cmap = 1;
opt.compress_cmap = false;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.show_counter = 0;
opt.delta_pair = 0.2*R;
opt.smallest_delta = 1e-3*R;
opt.ellipse_constant = true;
opt.use_tikhonov = true;
opt.tikhonov_tol = 1e-11;
opt.Nclust = 150;
opt.N_cmap = opt.N_c;
opt.use_interpolation = 'reduced_noconst';

panels = buildLaplaceCmapAlphaPanelPlan(opt);
n_panels = numel(panels);
n = 2*opt.N_cmap;
break_even_rank = floor(n*(sqrt(2)-1));

fallback_cache_file = fullfile(repo_root,'data', ...
    'sep28_capacitance_cmap_reference_free_rank_cache.mat');
results_file = fullfile(repo_root,'data', ...
    'sep28_capacitance_cmap_reference_free_rank_results.mat');
figure_file = fullfile(repo_root,'data', ...
    'sep28_capacitance_cmap_reference_free_rank.png');

fprintf(['Panels=%d, n=%d, basis q=%d, validation q=%d.\n' ...
    'Final tolerance %.1e; spatial budget %.1e; break-even rank %d.\n' ...
    'Core interpolation candidates: %s.\n\n'],n_panels,n,q_basis, ...
    q_validation,action_tolerance,spatial_action_tolerance, ...
    break_even_rank,mat2str(q_candidates));

%% Fixed canonical discretisation grids
N_c = opt.N_cmap;
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

%% Assemble every exact snapshot requested by the audit
to_delta = @(alpha) 2*R*(cosh(alpha)-1);
delta_chunks = cell(n_panels*(numel(q_candidates)+4),1);
ichunk = 0;
for ip = 1:n_panels
    panel = panels(ip);
    ichunk = ichunk+1;
    delta_chunks{ichunk} = to_delta(chebLobattoNodes( ...
        panel.alpha_lo,panel.alpha_hi,q_basis));
    ichunk = ichunk+1;
    delta_chunks{ichunk} = to_delta(chebLobattoNodes( ...
        panel.alpha_lo,panel.alpha_hi,q_validation));
    ichunk = ichunk+1;
    audit_alpha = panel.alpha_lo+(panel.alpha_hi-panel.alpha_lo)* ...
        audit_fraction(:);
    delta_chunks{ichunk} = to_delta(audit_alpha);
    ichunk = ichunk+1;
    delta_chunks{ichunk} = panel.reference_delta;
    for iq = 1:numel(q_candidates)
        ichunk = ichunk+1;
        delta_chunks{ichunk} = to_delta(chebLobattoNodes( ...
            panel.alpha_lo,panel.alpha_hi,q_candidates(iq)));
    end
end
delta_query = mergeNearlyEqual(vertcat(delta_chunks{1:ichunk}),1e-12);

fprintf('Requesting %d distinct exact Cmap snapshots...\n',numel(delta_query));
[C_pool,snapshot_source] = getExperimentSnapshots( ...
    delta_query,R,opt,grids,repo_root,fallback_cache_file);
fprintf(['  reused %d from production cache; requested %d from the ', ...
    'experiment fallback cache.\n\n'],snapshot_source.production_hits, ...
    snapshot_source.fallback_requests);

%% Compare centered and reference-free bases on each production panel
record_template = struct('panel',[],'base_index',[],'local_index',[], ...
    'delta_lo',[],'delta_hi',[],'production_q',[], ...
    'documented_centered_rank',[],'centered_spatial_rank',[], ...
    'centered_spatial_error',[],'direct_projection_rank_final',[], ...
    'direct_projection_error_final',[], ...
    'direct_projection_rank_spatial',[], ...
    'direct_projection_error_spatial',[], ...
    'full_interpolation_errors',[],'direct_end_to_end_ranks',[], ...
    'direct_end_to_end_errors',[],'selected_q',[], ...
    'selected_rank',[],'selected_error',[], ...
    'apply_cost',[],'apply_over_full',[],'speedup_over_full',[]);
records = repmat(record_template,n_panels,1);

for ip = 1:n_panels
    panel = panels(ip);
    basis_alpha = chebLobattoNodes( ...
        panel.alpha_lo,panel.alpha_hi,q_basis);
    validation_alpha = chebLobattoNodes( ...
        panel.alpha_lo,panel.alpha_hi,q_validation);
    audit_alpha = panel.alpha_lo+(panel.alpha_hi-panel.alpha_lo)* ...
        audit_fraction(:);

    C_basis = lookupMatrices(to_delta(basis_alpha),delta_query,C_pool);
    C_validation = lookupMatrices( ...
        to_delta(validation_alpha),delta_query,C_pool);
    C_audit = lookupMatrices(to_delta(audit_alpha),delta_query,C_pool);
    C_test = cat(3,C_validation,C_audit);
    test_norms = spectralNorms(C_test);
    C_ref = lookupMatrices(panel.reference_delta,delta_query,C_pool);
    C_ref = C_ref(:,:,1);

    direct_basis = commonBasisFromSnapshots(C_basis,n);
    centered_basis = commonBasisFromSnapshots(C_basis-C_ref,n);

    [direct_ranks,direct_errors] = selectRanksForTolerances( ...
        C_test,C_test,zeros(n),direct_basis.U,direct_basis.V, ...
        [action_tolerance spatial_action_tolerance], ...
        test_norms,rank_search_stride);
    [centered_rank,centered_error] = selectRanksForTolerances( ...
        C_test,C_test,C_ref,centered_basis.U,centered_basis.V, ...
        spatial_action_tolerance,test_norms,rank_search_stride);

    full_errors = nan(size(q_candidates));
    end_ranks = nan(size(q_candidates));
    end_errors = nan(size(q_candidates));
    for iq = 1:numel(q_candidates)
        q = q_candidates(iq);
        train_alpha = chebLobattoNodes(panel.alpha_lo,panel.alpha_hi,q);
        C_train = lookupMatrices(to_delta(train_alpha),delta_query,C_pool);
        held = heldOutMask(validation_alpha,train_alpha, ...
            panel.alpha_hi-panel.alpha_lo);
        C_interp_validation = evalMatrixChebBary( ...
            validation_alpha(held),train_alpha, ...
            chebBarycentricWeights(q),C_train);
        C_interp_audit = evalMatrixChebBary( ...
            audit_alpha,train_alpha,chebBarycentricWeights(q),C_train);
        C_interp = cat(3,C_interp_validation,C_interp_audit);
        C_exact = cat(3,C_validation(:,:,held),C_audit);
        exact_norms = spectralNorms(C_exact);
        full_errors(iq) = maxRelativeError( ...
            C_interp,C_exact,exact_norms);
        [end_ranks(iq),end_errors(iq)] = selectRanksForTolerances( ...
            C_interp,C_exact,zeros(n),direct_basis.U,direct_basis.V, ...
            action_tolerance,exact_norms,rank_search_stride);
    end

    selected_index = find(isfinite(end_ranks),1,'first');
    if isempty(selected_index)
        selected_q = NaN;
        selected_rank = NaN;
        selected_error = NaN;
        apply_cost = NaN;
        apply_over_full = NaN;
        speedup_over_full = NaN;
    else
        selected_q = q_candidates(selected_index);
        selected_rank = end_ranks(selected_index);
        selected_error = end_errors(selected_index);
        apply_cost = 2*n*selected_rank+selected_rank^2;
        apply_over_full = apply_cost/n^2;
        speedup_over_full = n^2/apply_cost;
    end

    records(ip).panel = ip;
    records(ip).base_index = panel.base_index;
    records(ip).local_index = panel.local_index;
    records(ip).delta_lo = panel.delta_lo;
    records(ip).delta_hi = panel.delta_hi;
    records(ip).production_q = panel.q;
    records(ip).documented_centered_rank = panel.spatial_rank;
    records(ip).centered_spatial_rank = centered_rank;
    records(ip).centered_spatial_error = centered_error;
    records(ip).direct_projection_rank_final = direct_ranks(1);
    records(ip).direct_projection_error_final = direct_errors(1);
    records(ip).direct_projection_rank_spatial = direct_ranks(2);
    records(ip).direct_projection_error_spatial = direct_errors(2);
    records(ip).full_interpolation_errors = full_errors;
    records(ip).direct_end_to_end_ranks = end_ranks;
    records(ip).direct_end_to_end_errors = end_errors;
    records(ip).selected_q = selected_q;
    records(ip).selected_rank = selected_rank;
    records(ip).selected_error = selected_error;
    records(ip).apply_cost = apply_cost;
    records(ip).apply_over_full = apply_over_full;
    records(ip).speedup_over_full = speedup_over_full;
end

%% Console report
fprintf('=== Reference-free panel ranks ===\n');
fprintf('%5s %5s %18s %5s %8s %8s %8s %9s %9s %9s\n', ...
    'panel','base','delta/R interval','qprod','r ctr','r C@1e-6', ...
    'r C@2.5e-7','q chosen','r end2end','cost/full');
for ip = 1:n_panels
    row = records(ip);
    interval = sprintf('[%.3g,%.3g]',row.delta_lo/R,row.delta_hi/R);
    fprintf('%5d %5d %18s %5d %8s %8s %8s %9s %9s %9s\n', ...
        row.panel,row.base_index,interval,row.production_q, ...
        numberOrDash(row.centered_spatial_rank,'%d'), ...
        numberOrDash(row.direct_projection_rank_final,'%d'), ...
        numberOrDash(row.direct_projection_rank_spatial,'%d'), ...
        numberOrDash(row.selected_q,'%d'), ...
        numberOrDash(row.selected_rank,'%d'), ...
        numberOrDash(row.apply_over_full,'%.3f'));
end

fprintf('\nEnd-to-end reference-free ranks by interpolation order:\n');
fprintf('%5s', 'panel');
for q = q_candidates
    fprintf('   q=%-3d',q);
end
fprintf('\n');
for ip = 1:n_panels
    fprintf('%5d',ip);
    for iq = 1:numel(q_candidates)
        fprintf(' %7s',numberOrDash( ...
            records(ip).direct_end_to_end_ranks(iq),'%d'));
    end
    fprintf('\n');
end

selected_ranks = [records.selected_rank];
finite_selected = isfinite(selected_ranks);
fprintf('\nBreak-even rank against one dense C apply: r <= %d.\n', ...
    break_even_rank);
if all(finite_selected)
    fprintf(['Selected end-to-end ranks: min %d, max %d; ', ...
        '%d/%d panels are below the break-even rank.\n'], ...
        min(selected_ranks),max(selected_ranks), ...
        sum(selected_ranks<=break_even_rank),n_panels);
    fprintf(['Reference-free apply/full cost ratio: %.3f--%.3f ', ...
        '(map-only speedup %.2f--%.2fx).\n'], ...
        min([records.apply_over_full]),max([records.apply_over_full]), ...
        min([records.speedup_over_full]),max([records.speedup_over_full]));
else
    fprintf('%d/%d panels found a passing q/r pair.\n', ...
        sum(finite_selected),n_panels);
end

%% Visible summary figure
f = figure('Name','Reference-free Cmap spatial ranks', ...
    'Color','w','Visible','on');
panel_index = 1:n_panels;

subplot(2,2,1);
plot(panel_index,[records.centered_spatial_rank],'-o','LineWidth',1.2, ...
    'DisplayName','centered, 2.5e-7');
hold on;
plot(panel_index,[records.direct_projection_rank_final],'-s', ...
    'LineWidth',1.2,'DisplayName','direct C, 1e-6');
plot(panel_index,[records.direct_projection_rank_spatial],'-d', ...
    'LineWidth',1.2,'DisplayName','direct C, 2.5e-7');
yline(break_even_rank,'--','DisplayName','dense-apply break-even');
xlabel('production panel');
ylabel('projection rank');
title('Exact-snapshot spatial projection');
legend('Location','best');
grid on;

subplot(2,2,2);
hold on;
for iq = 1:numel(q_candidates)
    rank_q = arrayfun(@(x) x.direct_end_to_end_ranks(iq),records);
    plot(panel_index,rank_q,'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('q=%d',q_candidates(iq)));
end
yline(break_even_rank,'--','DisplayName','dense-apply break-even');
xlabel('production panel');
ylabel('end-to-end rank');
title('After core interpolation');
legend('Location','best');
grid on;

subplot(2,2,3);
plot(panel_index,[records.apply_over_full],'-o','LineWidth',1.2);
hold on;
yline(1,'--','dense C cost');
xlabel('production panel');
ylabel('(2nr+r^2)/n^2');
title('Selected reference-free apply cost');
grid on;

subplot(2,2,4);
hold on;
for iq = 1:numel(q_candidates)
    errors_q = arrayfun(@(x) x.full_interpolation_errors(iq),records);
    semilogy(panel_index,errors_q,'-o','LineWidth',1.1, ...
        'DisplayName',sprintf('q=%d',q_candidates(iq)));
end
yline(action_tolerance,'--','DisplayName','target');
xlabel('production panel');
ylabel('relative action error');
title('Full interpolation floor before reduction');
legend('Location','best');
grid on;

saveas(f,figure_file);
drawnow;

save(results_file,'records','panels','q_basis','q_validation', ...
    'q_candidates','audit_fraction','action_tolerance', ...
    'spatial_action_tolerance','rank_search_stride','break_even_rank', ...
    'snapshot_source','opt','R','rng_seed');
fprintf('\nResults saved to %s\n',relativePath(results_file,repo_root));
fprintf('Figure saved to %s\n',relativePath(figure_file,repo_root));


%% Local helpers
function [C_query,source] = getExperimentSnapshots(delta_query,R,opt, ...
        grids,repo_root,fallback_cache_file)
exact_signature = buildLaplaceCmapInterpolationSignature(opt,false);
exact_signature = rmfield(exact_signature,'mode');
signature_id = laplaceInterpolationSignatureId(exact_signature);
production_file = fullfile(repo_root,'data', ...
    'laplace_cmap_interpolation',['exact_' signature_id '.mat']);

production_delta = zeros(0,1);
production_C = zeros(2*opt.N_cmap,2*opt.N_cmap,0);
if isfile(production_file)
    loaded = load(production_file,'cache');
    if isfield(loaded,'cache') && ...
            isequaln(loaded.cache.signature,exact_signature)
        production_delta = loaded.cache.delta;
        production_C = loaded.cache.C;
    end
end

nq = numel(delta_query);
C_query = zeros(2*opt.N_cmap,2*opt.N_cmap,nq);
found = false(nq,1);
for k = 1:nq
    index = findDelta(delta_query(k),production_delta);
    if ~isempty(index)
        C_query(:,:,k) = production_C(:,:,index);
        found(k) = true;
    end
end

missing_delta = delta_query(~found);
if ~isempty(missing_delta)
    [C_missing,~] = getOrBuildCmapSnapshots( ...
        missing_delta,R,opt,grids,fallback_cache_file);
    missing_index = find(~found);
    for k = 1:numel(missing_index)
        C_query(:,:,missing_index(k)) = C_missing(:,:,k);
    end
end

source = struct('production_file',production_file, ...
    'fallback_file',fallback_cache_file, ...
    'production_hits',sum(found), ...
    'fallback_requests',numel(missing_delta));
end

function index = findDelta(delta,pool)
index = find(abs(pool-delta) <= 1e-10*max(1,abs(delta)),1);
end

function [ranks,errors] = selectRanksForTolerances( ...
        C_source,C_exact,C_ref,U,V,tolerances,exact_norms,stride)
tolerances = tolerances(:).';
ranks = nan(size(tolerances));
errors = nan(size(tolerances));
r_max = min(size(U,2),size(V,2));
coarse_ranks = unique([1 stride:stride:r_max r_max]);
coarse_errors = nan(size(coarse_ranks));

for k = 1:numel(coarse_ranks)
    coarse_errors(k) = maxApproximationError( ...
        C_source,C_exact,C_ref,U,V,coarse_ranks(k),exact_norms);
end

for itol = 1:numel(tolerances)
    first_pass = find(coarse_errors <= tolerances(itol),1,'first');
    if isempty(first_pass)
        continue
    end
    if first_pass==1
        refine_ranks = coarse_ranks(1);
    else
        refine_ranks = (coarse_ranks(first_pass-1)+1): ...
            coarse_ranks(first_pass);
    end
    refine_errors = nan(size(refine_ranks));
    for k = 1:numel(refine_ranks)
        refine_errors(k) = maxApproximationError( ...
            C_source,C_exact,C_ref,U,V,refine_ranks(k),exact_norms);
    end
    local_pass = find(refine_errors <= tolerances(itol),1,'first');
    ranks(itol) = refine_ranks(local_pass);
    errors(itol) = refine_errors(local_pass);
end
end

function value = maxApproximationError( ...
        C_source,C_exact,C_ref,U,V,r,exact_norms)
U_r = U(:,1:r);
V_r = V(:,1:r);
value = 0;
for k = 1:size(C_exact,3)
    dC = C_source(:,:,k)-C_ref;
    approximation = C_ref+U_r*(U_r'*dC*V_r)*V_r';
    current = norm(approximation-C_exact(:,:,k),2)/exact_norms(k);
    value = max(value,current);
end
end

function values = spectralNorms(C)
values = zeros(size(C,3),1);
for k = 1:size(C,3)
    values(k) = max(norm(C(:,:,k),2),eps);
end
end

function value = maxRelativeError(C_approx,C_exact,exact_norms)
value = 0;
for k = 1:size(C_exact,3)
    value = max(value,norm(C_approx(:,:,k)-C_exact(:,:,k),2)/ ...
        exact_norms(k));
end
end

function mask = heldOutMask(validation_nodes,training_nodes,scale)
distance = min(abs(validation_nodes(:)-training_nodes(:)'),[],2);
mask = distance > 1e-10*max(1,abs(scale));
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

function p = relativePath(filename,root)
prefix = [root filesep];
if strncmp(filename,prefix,numel(prefix))
    p = filename(numel(prefix)+1:end);
else
    p = filename;
end
end
