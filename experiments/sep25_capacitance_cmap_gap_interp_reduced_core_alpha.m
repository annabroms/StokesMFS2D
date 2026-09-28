%% SEP25_CAPACITANCE_CMAP_GAP_INTERP_REDUCED_CORE_ALPHA
% Compare panelwise interpolation of the full Laplace capacitance pair map
% C(alpha) with interpolation of the smaller, panel-local core matrices
%
%   B(alpha) = U'*(C(alpha)-C_ref)*V,
%
% followed by C_hat(alpha)=C_ref+U*B_hat(alpha)*V'.  The comparison uses
% identical Chebyshev--Lobatto nodes for the direct and reduced methods and
% tests q=3, 5, 7, and 9 terms on each equal alpha panel.  Errors are the
% relative spectral/action errors of the final reconstructed C matrix, so
% the reduced result includes both spatial projection and B-interpolation
% error.  The core-only interpolation error is retained as a diagnostic.
%
% A q=33 grid supplies held-out validation points, and four independent
% random alpha values per panel provide a post-selection audit.  The basis
% rank is chosen using a tighter spatial error budget than the final target;
% otherwise projection alone can consume all of the 1e-6 allowance.
%
% Exact maps use Tikhonov regularisation, no Cmap compression, and a frozen
% ellipse discretisation with 150 nodes.  Figures are displayed directly
% and also saved under data/.
%
% Anna Broms, Sep 25, 2026

close all;
clearvars;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;
addpath(fullfile(repo_root,'experiments','cap_interp'));

fprintf('=== sep25_capacitance_cmap_gap_interp_reduced_core_alpha ===\n\n');

%% Configuration: reference_location is intentionally easy to change
R = 2;
P = 2;
interpolation_coordinate = 'alpha';
subdomain_counts = [1 2 4 8];
q_candidates = [3 5 7 9];
q_basis = 9;
q_validation = 33;
n_random_per_subdomain = 4;
action_tolerance = 1e-6;
spatial_action_tolerance = 2.5e-7;
reference_location = 'midpoint'; % 'minimum', 'midpoint', or 'maximum'
r_max = 160;
rank_search_stride = 4;
rng_seed = 20260925;

opt = getLaplace2Dparams(P,R);
opt.cmap = 1;
opt.compress_cmap = false;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.show_counter = 0;
opt.delta_pair = 0.2*R;
opt.ellipse_constant = true;
opt.use_tikhonov = true;
opt.tikhonov_tol = 1e-11;
opt.Nclust = 150;

delta_min = 1e-3*R;
delta_max = 0.999*opt.delta_pair;
delta_star = (R-opt.Rp_f)^2/opt.Rp_f;
opt.smallest_delta = delta_min;

to_alpha = @(delta) acosh(1+delta/(2*R));
to_delta = @(alpha) 2*R*(cosh(alpha)-1);
coordinate_name = 'alpha = acosh(1+delta/(2R))';

base_ranges = struct( ...
    'name',{'small-gap range','large-gap range'}, ...
    'delta_lo',{delta_min,delta_star}, ...
    'delta_hi',{delta_star,delta_max});

assert(any(strcmp(reference_location,{'minimum','midpoint','maximum'})), ...
    'Unsupported reference_location.');
assert(any(q_candidates==q_basis), ...
    'q_basis must be included in q_candidates.');

cache_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_subdomain_rank_cache.mat');
results_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_reduced_core_alpha_results.mat');
figure_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_reduced_core_alpha.png');

fprintf(['Coordinate: %s; frozen ellipse Nclust=%d at delta_min/R=%.3g.\n' ...
    'Direct C and reduced B use the same q candidates %s.\n' ...
    'Final target %.1e; spatial projection budget %.1e; reference %s.\n\n'], ...
    coordinate_name,opt.Nclust,opt.smallest_delta/R, ...
    mat2str(q_candidates),action_tolerance,spatial_action_tolerance, ...
    reference_location);

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

%% Plan equal alpha panels and the union of exact snapshots
rng(rng_seed);
n_plans = numel(base_ranges)*sum(subdomain_counts);
plan_template = struct('base_index',[],'n_subdomains',[], ...
    'local_index',[],'alpha_lo',[],'alpha_hi',[], ...
    'train_alpha',{{}},'basis_alpha',[],'validation_alpha',[], ...
    'random_alpha',[],'reference_delta',[]);
plans = repmat(plan_template,n_plans,1);
delta_chunks = cell(n_plans*(numel(q_candidates)+4),1);
ip = 0;
ichunk = 0;
for ib = 1:numel(base_ranges)
    alpha_base = to_alpha([base_ranges(ib).delta_lo; ...
        base_ranges(ib).delta_hi]);
    for ic = 1:numel(subdomain_counts)
        n_sub = subdomain_counts(ic);
        edges = linspace(alpha_base(1),alpha_base(2),n_sub+1);
        for isub = 1:n_sub
            ip = ip+1;
            plans(ip).base_index = ib;
            plans(ip).n_subdomains = n_sub;
            plans(ip).local_index = isub;
            plans(ip).alpha_lo = edges(isub);
            plans(ip).alpha_hi = edges(isub+1);
            plans(ip).train_alpha = cell(numel(q_candidates),1);
            for iq = 1:numel(q_candidates)
                q = q_candidates(iq);
                plans(ip).train_alpha{iq} = chebLobattoNodes( ...
                    plans(ip).alpha_lo,plans(ip).alpha_hi,q);
                ichunk = ichunk+1;
                delta_chunks{ichunk} = to_delta( ...
                    plans(ip).train_alpha{iq});
            end
            plans(ip).basis_alpha = chebLobattoNodes( ...
                plans(ip).alpha_lo,plans(ip).alpha_hi,q_basis);
            plans(ip).validation_alpha = chebLobattoNodes( ...
                plans(ip).alpha_lo,plans(ip).alpha_hi,q_validation);
            plans(ip).random_alpha = plans(ip).alpha_lo+ ...
                (plans(ip).alpha_hi-plans(ip).alpha_lo)* ...
                rand(n_random_per_subdomain,1);
            delta_lo = to_delta(plans(ip).alpha_lo);
            delta_hi = to_delta(plans(ip).alpha_hi);
            plans(ip).reference_delta = referenceDelta( ...
                delta_lo,delta_hi,reference_location);

            ichunk = ichunk+1;
            delta_chunks{ichunk} = to_delta(plans(ip).basis_alpha);
            ichunk = ichunk+1;
            delta_chunks{ichunk} = to_delta(plans(ip).validation_alpha);
            ichunk = ichunk+1;
            delta_chunks{ichunk} = to_delta(plans(ip).random_alpha);
            ichunk = ichunk+1;
            delta_chunks{ichunk} = plans(ip).reference_delta;
        end
    end
end
delta_pool = mergeNearlyEqual(vertcat(delta_chunks{1:ichunk}),1e-12);
fprintf('Building/reusing %d distinct exact Cmap snapshots...\n', ...
    numel(delta_pool));
[C_pool,~] = getOrBuildCmapSnapshots( ...
    delta_pool,R,opt,grids,cache_file);

%% Panel-local bases and interpolation-order comparison
local_template = struct('base_index',[],'n_subdomains',[], ...
    'local_index',[],'delta_lo',[],'delta_hi',[], ...
    'reference_delta',[],'spatial_rank',[], ...
    'projection_action_error',[],'direct_error_by_q',[], ...
    'reduced_error_by_q',[],'core_error_by_q',[], ...
    'direct_q_selected',[],'reduced_q_selected',[], ...
    'direct_held_error',[],'reduced_held_error',[], ...
    'direct_random_error',[],'reduced_random_error',[], ...
    'reduced_core_random_error',[]);
local = repmat(local_template,numel(plans),1);

for ip = 1:numel(plans)
    plan = plans(ip);
    C_ref_stack = lookupMatrices(plan.reference_delta,delta_pool,C_pool);
    C_ref = C_ref_stack(:,:,1);
    C_basis = lookupMatrices(to_delta(plan.basis_alpha),delta_pool,C_pool);
    C_validation = lookupMatrices( ...
        to_delta(plan.validation_alpha),delta_pool,C_pool);
    C_random = lookupMatrices(to_delta(plan.random_alpha),delta_pool,C_pool);

    basis = commonBasisFromSnapshots(C_basis-C_ref,r_max);
    held_basis = heldOutMask(plan.validation_alpha,plan.basis_alpha, ...
        plan.alpha_hi-plan.alpha_lo);
    C_rank_test = cat(3,C_validation(:,:,held_basis),C_random);
    dC_rank_test = C_rank_test-C_ref;
    [spatial_rank,projection_error] = selectSpatialRank( ...
        dC_rank_test,C_rank_test,basis.U,basis.V, ...
        spatial_action_tolerance,rank_search_stride);
    assert(isfinite(spatial_rank), ...
        'No spatial rank met the projection budget on panel %d.',ip);
    U = basis.U(:,1:spatial_rank);
    V = basis.V(:,1:spatial_rank);
    B_validation = projectCoreStack(C_validation-C_ref,U,V);
    B_random = projectCoreStack(C_random-C_ref,U,V);

    direct_error_by_q = nan(size(q_candidates));
    reduced_error_by_q = nan(size(q_candidates));
    core_error_by_q = nan(size(q_candidates));
    for iq = 1:numel(q_candidates)
        q = q_candidates(iq);
        alpha_train = plan.train_alpha{iq};
        C_train = lookupMatrices(to_delta(alpha_train),delta_pool,C_pool);
        held = heldOutMask(plan.validation_alpha,alpha_train, ...
            plan.alpha_hi-plan.alpha_lo);
        weights = chebBarycentricWeights(q);

        C_hat_direct = evalMatrixChebBary( ...
            plan.validation_alpha(held),alpha_train,weights,C_train);
        direct_error_by_q(iq) = max(relativeActionError( ...
            C_hat_direct,C_validation(:,:,held)));

        B_train = projectCoreStack(C_train-C_ref,U,V);
        B_hat = evalMatrixChebBary( ...
            plan.validation_alpha(held),alpha_train,weights,B_train);
        C_hat_reduced = reconstructFromCore(C_ref,U,V,B_hat);
        reduced_error_by_q(iq) = max(relativeActionError( ...
            C_hat_reduced,C_validation(:,:,held)));
        core_error_by_q(iq) = max(relativeCoreInterpolationError( ...
            B_hat,B_validation(:,:,held),C_validation(:,:,held)));
    end

    direct_index = find(direct_error_by_q <= action_tolerance,1,'first');
    reduced_index = find(reduced_error_by_q <= action_tolerance,1,'first');
    [direct_q,direct_held,direct_random] = selectedDirectAudit( ...
        direct_index,q_candidates,direct_error_by_q,plan, ...
        delta_pool,C_pool,C_random,to_delta);
    [reduced_q,reduced_held,reduced_random,reduced_core_random] = ...
        selectedReducedAudit(reduced_index,q_candidates, ...
        reduced_error_by_q,plan,delta_pool,C_pool,C_random,C_ref, ...
        U,V,B_random,to_delta);

    local(ip).base_index = plan.base_index;
    local(ip).n_subdomains = plan.n_subdomains;
    local(ip).local_index = plan.local_index;
    local(ip).delta_lo = to_delta(plan.alpha_lo);
    local(ip).delta_hi = to_delta(plan.alpha_hi);
    local(ip).reference_delta = plan.reference_delta;
    local(ip).spatial_rank = spatial_rank;
    local(ip).projection_action_error = projection_error;
    local(ip).direct_error_by_q = direct_error_by_q;
    local(ip).reduced_error_by_q = reduced_error_by_q;
    local(ip).core_error_by_q = core_error_by_q;
    local(ip).direct_q_selected = direct_q;
    local(ip).reduced_q_selected = reduced_q;
    local(ip).direct_held_error = direct_held;
    local(ip).reduced_held_error = reduced_held;
    local(ip).direct_random_error = direct_random;
    local(ip).reduced_random_error = reduced_random;
    local(ip).reduced_core_random_error = reduced_core_random;
end

%% Aggregate all panels in each base range
summary_template = struct('base_index',[],'base_name','', ...
    'n_subdomains',[],'spatial_ranks',[],'direct_q_local',[], ...
    'reduced_q_local',[],'direct_q_max',[],'reduced_q_max',[], ...
    'direct_held_error',[],'reduced_held_error',[], ...
    'direct_random_error',[],'reduced_random_error',[], ...
    'reduced_core_random_error',[],'max_projection_error',[], ...
    'direct_validation_pass',[],'reduced_validation_pass',[], ...
    'direct_random_pass',[],'reduced_random_pass',[]);
summary = repmat(summary_template, ...
    numel(base_ranges)*numel(subdomain_counts),1);
isummary = 0;
for ib = 1:numel(base_ranges)
    for ic = 1:numel(subdomain_counts)
        n_sub = subdomain_counts(ic);
        take = find([local.base_index]==ib & ...
            [local.n_subdomains]==n_sub);
        [~,order] = sort([local(take).local_index]);
        take = take(order);
        direct_q = [local(take).direct_q_selected];
        reduced_q = [local(take).reduced_q_selected];
        direct_ok = all(isfinite(direct_q));
        reduced_ok = all(isfinite(reduced_q));

        isummary = isummary+1;
        summary(isummary).base_index = ib;
        summary(isummary).base_name = base_ranges(ib).name;
        summary(isummary).n_subdomains = n_sub;
        summary(isummary).spatial_ranks = [local(take).spatial_rank];
        summary(isummary).direct_q_local = direct_q;
        summary(isummary).reduced_q_local = reduced_q;
        summary(isummary).direct_q_max = finiteMaximum(direct_q,direct_ok);
        summary(isummary).reduced_q_max = finiteMaximum(reduced_q,reduced_ok);
        summary(isummary).direct_held_error = finiteMaximum( ...
            [local(take).direct_held_error],direct_ok);
        summary(isummary).reduced_held_error = finiteMaximum( ...
            [local(take).reduced_held_error],reduced_ok);
        summary(isummary).direct_random_error = finiteMaximum( ...
            [local(take).direct_random_error],direct_ok);
        summary(isummary).reduced_random_error = finiteMaximum( ...
            [local(take).reduced_random_error],reduced_ok);
        summary(isummary).reduced_core_random_error = finiteMaximum( ...
            [local(take).reduced_core_random_error],reduced_ok);
        summary(isummary).max_projection_error = max( ...
            [local(take).projection_action_error]);
        summary(isummary).direct_validation_pass = direct_ok;
        summary(isummary).reduced_validation_pass = reduced_ok;
        summary(isummary).direct_random_pass = direct_ok && ...
            summary(isummary).direct_random_error <= action_tolerance;
        summary(isummary).reduced_random_pass = reduced_ok && ...
            summary(isummary).reduced_random_error <= action_tolerance;
    end
end

%% Report
fprintf('\n=== Full C versus reduced B interpolation at %.1e ===\n', ...
    action_tolerance);
fprintf(['%-17s %5s %-19s %-19s %-17s %10s %10s ', ...
    '%10s %10s\n'], ...
    'base range','parts','spatial ranks','q: direct C','q: reduced B', ...
    'C held','C random','B held','B random');
for k = 1:numel(summary)
    row = summary(k);
    fprintf(['%-17s %5d %-19s %-19s %-17s %10s %10s ', ...
        '%10s %10s\n'],row.base_name,row.n_subdomains, ...
        mat2str(row.spatial_ranks),mat2str(row.direct_q_local), ...
        mat2str(row.reduced_q_local), ...
        numberOrDash(row.direct_held_error,'%.3e'), ...
        numberOrDash(row.direct_random_error,'%.3e'), ...
        numberOrDash(row.reduced_held_error,'%.3e'), ...
        numberOrDash(row.reduced_random_error,'%.3e'));
end

fprintf(['\nB held/random are end-to-end C errors after reconstruction; ', ...
    'they include the panel-local spatial projection error.\n']);
fprintf(['The smaller core does not by itself lower polynomial degree: ', ...
    'fixed linear projection and barycentric interpolation commute.\n']);
fprintf(['Any reduction in selected q therefore comes from discarding ', ...
    'directions below the spatial error budget.\n\n']);

%% Visible figure
f = figure('Name','Direct C versus reduced B interpolation in alpha', ...
    'Color','w','Visible','on');
for ib = 1:numel(base_ranges)
    take = find([summary.base_index]==ib);
    counts = [summary(take).n_subdomains];
    subplot(2,numel(base_ranges),ib);
    plot(counts,[summary(take).direct_q_max],'-o','LineWidth',1.2, ...
        'DisplayName','direct C');
    hold on;
    plot(counts,[summary(take).reduced_q_max],'-s','LineWidth',1.2, ...
        'DisplayName','reduced B');
    xlabel('equal alpha panels');
    ylabel('maximum local q');
    title(base_ranges(ib).name);
    xticks(subdomain_counts);
    ylim([2 10]);
    legend('Location','best');
    grid on;

    subplot(2,numel(base_ranges),ib+numel(base_ranges));
    semilogy(counts,[summary(take).direct_random_error],'-o', ...
        'LineWidth',1.2,'DisplayName','direct C random');
    hold on;
    semilogy(counts,[summary(take).reduced_random_error],'-s', ...
        'LineWidth',1.2,'DisplayName','reduced B random');
    yline(action_tolerance,'--','DisplayName','target');
    xlabel('equal alpha panels');
    ylabel('relative action error');
    title('post-selection audit');
    xticks(subdomain_counts);
    legend('Location','best');
    grid on;
end
saveas(f,figure_file);
drawnow;

save(results_file,'local','summary','plans','base_ranges', ...
    'subdomain_counts','q_candidates','q_basis','q_validation', ...
    'n_random_per_subdomain','action_tolerance', ...
    'spatial_action_tolerance','reference_location','r_max', ...
    'rank_search_stride','rng_seed','interpolation_coordinate', ...
    'coordinate_name','delta_star','opt','R','-v7.3');
fprintf('Results saved to %s\n',relativePath(results_file,repo_root));
fprintf('Figure saved to %s\n',relativePath(figure_file,repo_root));


%% Local helpers
function delta_ref = referenceDelta(delta_lo,delta_hi,location)
switch location
    case 'minimum'
        delta_ref = delta_lo;
    case 'midpoint'
        delta_ref = (delta_lo+delta_hi)/2;
    case 'maximum'
        delta_ref = delta_hi;
    otherwise
        error('Unknown reference location "%s".',location);
end
end

function B = projectCoreStack(dC,U,V)
m = size(dC,3);
r = size(U,2);
B = zeros(r,r,m);
for k = 1:m
    B(:,:,k) = U'*dC(:,:,k)*V;
end
end

function C = reconstructFromCore(C_ref,U,V,B)
n = size(C_ref,1);
m = size(B,3);
C = zeros(n,n,m);
for k = 1:m
    C(:,:,k) = C_ref+U*B(:,:,k)*V';
end
end

function error = relativeActionError(C_hat,C_exact)
m = size(C_exact,3);
error = zeros(m,1);
for k = 1:m
    error(k) = norm(C_hat(:,:,k)-C_exact(:,:,k),2)/ ...
        max(norm(C_exact(:,:,k),2),eps);
end
end

function error = relativeCoreInterpolationError(B_hat,B_exact,C_exact)
m = size(B_exact,3);
error = zeros(m,1);
for k = 1:m
    error(k) = norm(B_hat(:,:,k)-B_exact(:,:,k),2)/ ...
        max(norm(C_exact(:,:,k),2),eps);
end
end

function [q_selected,held_error,random_error] = selectedDirectAudit( ...
        index,q_candidates,error_by_q,plan,delta_pool,C_pool,C_random, ...
        to_delta)
if isempty(index)
    q_selected = NaN;
    held_error = NaN;
    random_error = NaN;
    return
end
q_selected = q_candidates(index);
held_error = error_by_q(index);
alpha_train = plan.train_alpha{index};
C_train = lookupMatrices(to_delta(alpha_train),delta_pool,C_pool);
C_hat = evalMatrixChebBary(plan.random_alpha,alpha_train, ...
    chebBarycentricWeights(q_selected),C_train);
random_error = max(relativeActionError(C_hat,C_random));
end

function [q_selected,held_error,random_error,core_random_error] = ...
        selectedReducedAudit(index,q_candidates,error_by_q,plan, ...
        delta_pool,C_pool,C_random,C_ref,U,V,B_random,to_delta)
if isempty(index)
    q_selected = NaN;
    held_error = NaN;
    random_error = NaN;
    core_random_error = NaN;
    return
end
q_selected = q_candidates(index);
held_error = error_by_q(index);
alpha_train = plan.train_alpha{index};
C_train = lookupMatrices(to_delta(alpha_train),delta_pool,C_pool);
B_train = projectCoreStack(C_train-C_ref,U,V);
B_hat = evalMatrixChebBary(plan.random_alpha,alpha_train, ...
    chebBarycentricWeights(q_selected),B_train);
C_hat = reconstructFromCore(C_ref,U,V,B_hat);
random_error = max(relativeActionError(C_hat,C_random));
core_random_error = max(relativeCoreInterpolationError( ...
    B_hat,B_random,C_random));
end

function [r_selected,selected_error] = selectSpatialRank( ...
        dC_test,C_test,U,V,tolerance,stride)
r_max_local = size(U,2);
coarse_ranks = unique([1,stride:stride:r_max_local,r_max_local]);
coarse_error = nan(size(coarse_ranks));
first_pass = [];
for k = 1:numel(coarse_ranks)
    coarse_error(k) = maxProjectionActionError( ...
        dC_test,C_test,U,V,coarse_ranks(k));
    if coarse_error(k) <= tolerance
        first_pass = k;
        break
    end
end
if isempty(first_pass)
    r_selected = NaN;
    selected_error = NaN;
    return
end
if first_pass==1
    refine_ranks = coarse_ranks(1);
else
    refine_ranks = (coarse_ranks(first_pass-1)+1): ...
        coarse_ranks(first_pass);
end
refine_error = nan(size(refine_ranks));
for k = 1:numel(refine_ranks)
    refine_error(k) = maxProjectionActionError( ...
        dC_test,C_test,U,V,refine_ranks(k));
end
first_refined_pass = find(refine_error <= tolerance,1,'first');
r_selected = refine_ranks(first_refined_pass);
selected_error = refine_error(first_refined_pass);
end

function error = maxProjectionActionError(dC,C,U,V,r)
U_r = U(:,1:r);
V_r = V(:,1:r);
error = 0;
for k = 1:size(C,3)
    projected = U_r*(U_r'*dC(:,:,k)*V_r)*V_r';
    current = norm(dC(:,:,k)-projected,2)/ ...
        max(norm(C(:,:,k),2),eps);
    error = max(error,current);
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

function value = finiteMaximum(values,valid)
if valid
    value = max(values);
else
    value = NaN;
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
