%% SEP25_CAPACITANCE_CMAP_GAP_INTERP_SUBDOMAINS Piecewise interpolation
% of the Tikhonov-regularised Laplace capacitance pair map in the gap.
%
% This experiment extends sep25_capacitance_cmap_gap_interp.m.  The gap
% range is divided at the former ellipse-activation distance.  Each
% of these two BASE RANGES is additionally split into 1, 2, 4, or 8 equal
% numerical subdomains in the selected interpolation coordinate.  The
% default coordinate is s=log(delta/R).  The companion
% sep25_capacitance_cmap_gap_interp_subdomains_alpha.m wrapper repeats the
% same experiment with alpha=acosh(1+delta/(2R)).
%
% For each local interpolant, q training nodes are tested against
% exact matrices at the held-out nodes of a q=33 nested grid.  The script
% selects the smallest q whose relative spectral/action error
%
%     ||C(delta)-Chat(delta)||_2 / ||C(delta)||_2
%
% is below a prescribed tolerance.  This is an a posteriori validation
% criterion, not an internal tolerance of the barycentric formula.  Four
% additional random gaps per numerical subdomain provide an independent
% check after q has been selected.
%
% The primary quantity is q_max, the largest number of terms in any local
% interpolation sum.  The total number n_unique of distinct exact training
% snapshots is reported only as secondary offline bookkeeping; it is not an
% optimization criterion.  Shared endpoints are counted once in n_unique,
% although each local barycentric formula still contains its own endpoint
% term.
%
% Tikhonov regularisation is used for the fine-pair and peanut
% pseudoinverses.  The ellipse construction is frozen at delta_min with
% ellipse_constant=true and Nclust=150, so its nodes do not enter or leave
% as the gap changes.  Cmap compression is disabled because compress_cmap
% is a separate TSVD operation.  Cached exact matrices and figures are
% written under data/, which is git-ignored.
%
% Anna Broms, Sep 25, 2026

if exist('cap_interp_coordinate_override','var') && ...
        ~isempty(cap_interp_coordinate_override)
    interpolation_coordinate = cap_interp_coordinate_override;
else
    interpolation_coordinate = 'log_delta';
end
close all;
clearvars -except interpolation_coordinate;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;
addpath(fullfile(repo_root,'experiments','cap_interp'));

fprintf('=== sep25_capacitance_cmap_gap_interp_subdomains ===\n\n');

%% Configuration
R = 2;
P = 2;

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

% Equal numerical subdivisions in the selected coordinate inside each
% base gap range.  The number of offline snapshots is not minimised; the
% objective is a short local barycentric sum.
subdomain_counts = [1 2 4 8];
q_candidates = [3 5 9];
q_validation = 33;
action_tolerances = [1e-4 1e-5 3e-6 1e-6];
primary_action_tol = 1e-6;
n_random_per_subdomain = 4;
rng_seed = 20260925;

assert(any(action_tolerances==primary_action_tol), ...
    'primary_action_tol must be included in action_tolerances.');
nesting_ratios = (q_validation-1)./(q_candidates-1);
assert(all(abs(nesting_ratios-round(nesting_ratios)) < 10*eps), ...
    'Every candidate grid must be nested in the q_validation grid.');

switch interpolation_coordinate
    case 'log_delta'
        coordinate_name = 's = log(delta/R)';
        coordinate_tag = 'log_delta';
        to_param = @(delta) log(delta/R);
        to_delta = @(s) R*exp(s);
    case 'alpha'
        coordinate_name = 'alpha = acosh(1+delta/(2R))';
        coordinate_tag = 'alpha';
        to_param = @(delta) acosh(1+delta/(2*R));
        to_delta = @(alpha) 2*R*(cosh(alpha)-1);
    otherwise
        error('Unknown interpolation coordinate "%s".', ...
            interpolation_coordinate);
end

cache_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_tikhonov_constant150_cache.mat');
if strcmp(interpolation_coordinate,'log_delta')
    results_name = 'sep25_capacitance_cmap_gap_interp_subdomains_results.mat';
    figure_name = 'sep25_capacitance_cmap_gap_interp_subdomains.png';
else
    results_name = sprintf(['sep25_capacitance_cmap_gap_interp_', ...
        'subdomains_%s_results.mat'],coordinate_tag);
    figure_name = sprintf(['sep25_capacitance_cmap_gap_interp_', ...
        'subdomains_%s.png'],coordinate_tag);
end
results_file = fullfile(repo_root,'data',results_name);
figure_file = fullfile(repo_root,'data',figure_name);

fprintf(['Coordinate: %s.  Frozen ellipse: Nclust=%d at delta_min/R=%.3g.\n' ...
    'Tikhonov eta=lambda/sigma_max=%.1e; primary action tolerance=%.1e.\n' ...
    'Testing subdomain counts %s, q candidates %s, q_validation=%d.\n\n'], ...
    coordinate_name,opt.Nclust,opt.smallest_delta/R, ...
    opt.tikhonov_tol,primary_action_tol,mat2str(subdomain_counts), ...
    mat2str(q_candidates),q_validation);

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

%% Two base gap ranges; frozen ellipse sources are active in both
% delta_star is retained as a reproducible comparison split with the older
% gap-following discretisation.  With ellipse_constant=true it is not an
% activation boundary and no transition bracket is required.
physical_regimes = struct( ...
    'name',{'small-gap range','large-gap range'}, ...
    'delta_lo',{delta_min,delta_star}, ...
    'delta_hi',{delta_star,delta_max});

%% Plan all exact samples before building them
% Building the union in one call avoids repeatedly saving an expanding
% cache file.  Candidate training grids are nested subsets of each q=33
% validation grid, so no separate construction is needed for them.
rng(rng_seed);
plans = struct([]);
all_delta = zeros(0,1);
for ip = 1:numel(physical_regimes)
    regime = physical_regimes(ip);
    s_regime = [to_param(regime.delta_lo),to_param(regime.delta_hi)];
    for ic = 1:numel(subdomain_counts)
        n_sub = subdomain_counts(ic);
        edges = linspace(s_regime(1),s_regime(2),n_sub+1);
        for isub = 1:n_sub
            plan.physical_index = ip;
            plan.n_subdomains = n_sub;
            plan.local_index = isub;
            plan.s_lo = edges(isub);
            plan.s_hi = edges(isub+1);
            plan.validation_s = chebLobattoNodes( ...
                plan.s_lo,plan.s_hi,q_validation);
            plan.random_s = plan.s_lo+(plan.s_hi-plan.s_lo)* ...
                rand(n_random_per_subdomain,1);
            plans = [plans; plan]; %#ok<AGROW>
            all_delta = [all_delta; to_delta(plan.validation_s); ...
                to_delta(plan.random_s)]; %#ok<AGROW>
        end
    end
end

delta_pool = mergeNearlyEqual(all_delta,1e-12);
fprintf('Building/reusing %d distinct exact Cmap snapshots...\n',numel(delta_pool));
[C_pool,~] = getOrBuildCmapSnapshots(delta_pool,R,opt,grids,cache_file);

%% Local interpolation-order selection
local = repmat(struct('physical_index',[],'n_subdomains',[], ...
    'local_index',[],'s_lo',[],'s_hi',[],'error_by_q',[], ...
    'q_selected',[],'selected_validation_error',[], ...
    'random_error',[]),numel(plans),1);

for k = 1:numel(plans)
    plan = plans(k);
    C_validation = lookupMatrices(to_delta(plan.validation_s),delta_pool,C_pool);
    error_by_q = nan(size(q_candidates));

    for iq = 1:numel(q_candidates)
        q = q_candidates(iq);
        s_train = chebLobattoNodes(plan.s_lo,plan.s_hi,q);
        C_train = lookupMatrices(to_delta(s_train),delta_pool,C_pool);
        is_held = heldOutMask(plan.validation_s,s_train,plan.s_hi-plan.s_lo);
        C_hat = evalMatrixChebBary(plan.validation_s(is_held),s_train, ...
            chebBarycentricWeights(q),C_train);
        error_by_q(iq) = max(relativeActionError( ...
            C_hat,C_validation(:,:,is_held)));
    end

    q_selected = nan(size(action_tolerances));
    selected_error = nan(size(action_tolerances));
    for it = 1:numel(action_tolerances)
        iq = find(error_by_q <= action_tolerances(it),1,'first');
        if ~isempty(iq)
            q_selected(it) = q_candidates(iq);
            selected_error(it) = error_by_q(iq);
        end
    end

    iprimary = find(action_tolerances==primary_action_tol,1);
    q_primary = q_selected(iprimary);
    if isfinite(q_primary)
        s_train = chebLobattoNodes(plan.s_lo,plan.s_hi,q_primary);
        C_train = lookupMatrices(to_delta(s_train),delta_pool,C_pool);
        C_random = lookupMatrices(to_delta(plan.random_s),delta_pool,C_pool);
        C_hat_random = evalMatrixChebBary(plan.random_s,s_train, ...
            chebBarycentricWeights(q_primary),C_train);
        random_error = max(relativeActionError(C_hat_random,C_random));
    else
        random_error = NaN;
    end

    local(k).physical_index = plan.physical_index;
    local(k).n_subdomains = plan.n_subdomains;
    local(k).local_index = plan.local_index;
    local(k).s_lo = plan.s_lo;
    local(k).s_hi = plan.s_hi;
    local(k).error_by_q = error_by_q;
    local(k).q_selected = q_selected;
    local(k).selected_validation_error = selected_error;
    local(k).random_error = random_error;
end

%% Aggregate local results, counting shared endpoints only once
summary = struct([]);
for ip = 1:numel(physical_regimes)
    for ic = 1:numel(subdomain_counts)
        n_sub = subdomain_counts(ic);
        take = find([local.physical_index]==ip & ...
            [local.n_subdomains]==n_sub);
        for it = 1:numel(action_tolerances)
            q_local = arrayfun(@(x) x.q_selected(it),local(take));
            ok = all(isfinite(q_local));
            if ok
                train_s = zeros(0,1);
                for j = 1:numel(take)
                    item = local(take(j));
                    train_s = [train_s; chebLobattoNodes( ...
                        item.s_lo,item.s_hi,q_local(j))]; %#ok<AGROW>
                end
                n_unique = numel(mergeNearlyEqual(train_s,1e-12));
                n_terms_local_sum = sum(q_local);
                q_max = max(q_local);
                validation_error = max(arrayfun(@(x) ...
                    x.selected_validation_error(it),local(take)));
            else
                n_unique = NaN;
                n_terms_local_sum = NaN;
                q_max = NaN;
                validation_error = NaN;
            end

            row.physical_index = ip;
            row.physical_name = physical_regimes(ip).name;
            row.n_subdomains = n_sub;
            row.action_tolerance = action_tolerances(it);
            row.pass = ok;
            row.q_local = q_local;
            row.q_max = q_max;
            row.n_terms_local_sum = n_terms_local_sum;
            row.n_unique_training_nodes = n_unique;
            row.validation_error = validation_error;
            if action_tolerances(it)==primary_action_tol && ok
                row.random_error = max([local(take).random_error]);
            else
                row.random_error = NaN;
            end
            summary = [summary; row]; %#ok<AGROW>
        end
    end
end

%% Report
fprintf('\n=== Smallest local q selected from %s ===\n',mat2str(q_candidates));
fprintf('%-18s %5s %9s %-18s %6s %8s %11s %11s\n', ...
    'base gap range','pieces','tolerance','q on pieces','q_max', ...
    'unique','held-out','random');
for k = 1:numel(summary)
    r = summary(k);
    fprintf('%-18s %5d %9.1e %-18s %6s %8s %11s %11s\n', ...
        r.physical_name,r.n_subdomains,r.action_tolerance, ...
        mat2str(r.q_local),numberOrDash(r.q_max,'%d'), ...
        numberOrDash(r.n_unique_training_nodes,'%d'), ...
        numberOrDash(r.validation_error,'%.3e'), ...
        numberOrDash(r.random_error,'%.3e'));
end

fprintf(['\nq is the number of barycentric terms (polynomial degree q-1) ', ...
    'on one numerical subdomain.\n']);
fprintf(['unique is the total number of distinct exact training snapshots ', ...
    'for the whole physical regime.\n']);
fprintf(['Certification uses the non-training nodes of a nested q=%d grid; ', ...
    'these validation matrices are not counted as online interpolant ', ...
    'terms.\n\n'],q_validation);

printPrimaryComparison(summary,physical_regimes,subdomain_counts, ...
    primary_action_tol);

%% Plot the primary-tolerance tradeoff
f = figure('Name',['Cmap subdivision tradeoff: ' coordinate_name], ...
    'Color','w','Visible','on');
for ip = 1:numel(physical_regimes)
    take = find([summary.physical_index]==ip & ...
        [summary.action_tolerance]==primary_action_tol);
    subplot(2,2,ip);
    plot([summary(take).n_subdomains],[summary(take).q_max],'-o', ...
        'LineWidth',1.2);
    xlabel('numerical subdomains'); ylabel('maximum local q');
    title([physical_regimes(ip).name ': local order']); grid on;

    subplot(2,2,ip+2);
    semilogy([summary(take).n_subdomains], ...
        [summary(take).validation_error],'-o','LineWidth',1.2, ...
        'DisplayName','held-out'); hold on;
    semilogy([summary(take).n_subdomains], ...
        [summary(take).random_error],'-s','LineWidth',1.2, ...
        'DisplayName','random');
    yline(primary_action_tol,'--','DisplayName','target');
    xlabel('numerical subdomains'); ylabel('relative action error');
    title([physical_regimes(ip).name ': accuracy']);
    legend('Location','best'); grid on;
end
saveas(f,figure_file);
drawnow;

save(results_file,'summary','local','plans','physical_regimes', ...
    'subdomain_counts','q_candidates','q_validation','action_tolerances', ...
    'primary_action_tol','n_random_per_subdomain','rng_seed','delta_star', ...
    'interpolation_coordinate','coordinate_name','opt','R','-v7.3');
fprintf('Results saved to %s\n',relativePath(results_file,repo_root));
fprintf('Figure saved to %s\n',relativePath(figure_file,repo_root));


%% Local helpers
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

function printPrimaryComparison(summary,physical_regimes,counts,tolerance)
fprintf('=== Primary tolerance %.1e ===\n',tolerance);
for ip = 1:numel(physical_regimes)
    take = find([summary.physical_index]==ip & ...
        [summary.action_tolerance]==tolerance);
    one = summary(take([summary(take).n_subdomains]==1));
    fprintf('%s:\n',physical_regimes(ip).name);
    for ic = 1:numel(counts)
        r = summary(take([summary(take).n_subdomains]==counts(ic)));
        if r.pass
            fprintf(['  %d piece(s): q=%s, q_max=%d, %d distinct training ', ...
                'nodes, held-out %.3e, random %.3e'], ...
                counts(ic),mat2str(r.q_local),r.q_max, ...
                r.n_unique_training_nodes,r.validation_error,r.random_error);
            if one.pass
                fprintf(' (local-q change versus one piece: %+d)', ...
                    r.q_max-one.q_max);
            end
            fprintf('\n');
        else
            fprintf('  %d piece(s): no q in the tested set passed\n',counts(ic));
        end
    end
end
fprintf('\n');
end

function p = relativePath(filename,root)
prefix = [root filesep];
if strncmp(filename,prefix,numel(prefix))
    p = filename(numel(prefix)+1:end);
else
    p = filename;
end
end
