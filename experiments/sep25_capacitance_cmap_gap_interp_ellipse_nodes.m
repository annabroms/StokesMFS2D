%% SEP25_CAPACITANCE_CMAP_GAP_INTERP_ELLIPSE_NODES Ellipse-density sweep
% for piecewise interpolation of the Tikhonov-regularised capacitance Cmap.
%
% The eight-piece experiment in
% sep25_capacitance_cmap_gap_interp_subdomains.m reaches a 1e-5 action-error
% target with five barycentric terms on every numerical subdomain, but the
% ellipse-enhanced physical regime does not pass a 1e-6 target everywhere.
% This focused experiment tests whether a modest increase in the ellipse
% discretisation parameter opt.Nclust removes that error floor.
%
% Only the ellipse-enhanced physical interval is tested.  It is divided
% uniformly into eight numerical subdomains in s=log(delta/R).  For each
% Nclust in [100 110 120 130], local q=5,9,17 interpolants are validated at
% the nontraining points of a nested q=33 exact grid.  The smallest passing
% q is selected independently on every subdomain, followed by four random
% validation gaps per subdomain.  Figures are shown directly and a copy is
% also saved under data/.
%
% The error criterion is
%
%   ||C(delta)-Chat(delta)||_2 / ||C(delta)||_2 <= 1e-6.
%
% Anna Broms, Sep 25, 2026

close all;
clearvars;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;
addpath(fullfile(repo_root,'experiments','cap_interp'));

fprintf('=== sep25_capacitance_cmap_gap_interp_ellipse_nodes ===\n\n');

%% Configuration
R = 2;
P = 2;
opt_base = getLaplace2Dparams(P,R);
opt_base.cmap = 1;
opt_base.compress_cmap = false;
opt_base.reuse_pair_basis_by_sep = false;
opt_base.parallel_precomp = false;
opt_base.show_counter = 0;
opt_base.delta_pair = 0.2*R;
opt_base.ellipse_constant = false;
opt_base.use_tikhonov = true;
opt_base.tikhonov_tol = 1e-11;

Nclust_values = [100 110 120 130];
n_subdomains = 8;
q_candidates = [5 9 17];
q_validation = 33;
action_tolerance = 1e-6;
n_random_per_subdomain = 4;
rng_seed = 20260925;
transition_margin = 1e-8;

delta_min = 1e-3*R;
delta_star = (R-opt_base.Rp_f)^2/opt_base.Rp_f;
delta_hi = delta_star*(1-transition_margin);
to_param = @(delta) log(delta/R);
to_delta = @(s) R*exp(s);

cache_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_tikhonov_cache.mat');
results_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_ellipse_nodes_results.mat');
figure_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_ellipse_nodes.png');

fprintf(['Target %.1e, Nclust=%s, %d numerical subdomains, ', ...
    'q candidates=%s.\n\n'],action_tolerance,mat2str(Nclust_values), ...
    n_subdomains,mat2str(q_candidates));

%% Fixed discretisation grids
N_c = opt_base.N_c;
N_f = opt_base.N_f;
nout_c = ceil(opt_base.a_c*N_c);
tout_c = linspace(0,2*pi,nout_c+1)';
tout_c = tout_c(1:end-1);
grids.rbase_out_c = R*(cos(tout_c)+1i*sin(tout_c));

tin_c = linspace(0,2*pi,N_c+1)';
tin_c = tin_c(1:end-1);
grids.rbase_in_c = opt_base.Rp_c*(cos(tin_c)+1i*sin(tin_c));

tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
grids.rbase_in_f = opt_base.Rp_f*(cos(tin_f)+1i*sin(tin_f));

nout_f = ceil(opt_base.a_f*N_f);
tout_f = linspace(0,2*pi,nout_f+1)';
tout_f = tout_f(1:end-1);
grids.rout_base_f = R*(cos(tout_f)+1i*sin(tout_f));

%% Use the same validation and random gaps for every ellipse density
s_edges = linspace(to_param(delta_min),to_param(delta_hi),n_subdomains+1);
rng(rng_seed);
plans = repmat(struct('s_lo',[],'s_hi',[],'validation_s',[], ...
    'random_s',[]),n_subdomains,1);
all_delta = zeros(0,1);
for isub = 1:n_subdomains
    plans(isub).s_lo = s_edges(isub);
    plans(isub).s_hi = s_edges(isub+1);
    plans(isub).validation_s = chebLobattoNodes( ...
        plans(isub).s_lo,plans(isub).s_hi,q_validation);
    plans(isub).random_s = plans(isub).s_lo+ ...
        (plans(isub).s_hi-plans(isub).s_lo)* ...
        rand(n_random_per_subdomain,1);
    all_delta = [all_delta; to_delta(plans(isub).validation_s); ...
        to_delta(plans(isub).random_s)]; %#ok<AGROW>
end
delta_pool = mergeNearlyEqual(all_delta,1e-12);

%% Sweep opt.Nclust
results = repmat(struct('Nclust',[],'pass',[],'q_local',[],'q_max',[], ...
    'held_out_error',[],'random_error',[],'qmax_test_error',[], ...
    'retained_ellipse_nodes',[],'ellipse_count_changes',[], ...
    'local_error_by_q',[]),numel(Nclust_values),1);

for in = 1:numel(Nclust_values)
    opt = opt_base;
    opt.Nclust = Nclust_values(in);
    fprintf('--- Nclust=%d: building/reusing %d exact snapshots ---\n', ...
        opt.Nclust,numel(delta_pool));
    [C_pool,info_pool] = getOrBuildCmapSnapshots( ...
        delta_pool,R,opt,grids,cache_file);

    q_local = nan(1,n_subdomains);
    held_local = nan(1,n_subdomains);
    random_local = nan(1,n_subdomains);
    qmax_local = nan(1,n_subdomains);
    error_by_q = nan(n_subdomains,numel(q_candidates));

    for isub = 1:n_subdomains
        plan = plans(isub);
        C_validation = lookupMatrices( ...
            to_delta(plan.validation_s),delta_pool,C_pool);

        for iq = 1:numel(q_candidates)
            q = q_candidates(iq);
            s_train = chebLobattoNodes(plan.s_lo,plan.s_hi,q);
            C_train = lookupMatrices(to_delta(s_train),delta_pool,C_pool);
            held = heldOutMask(plan.validation_s,s_train, ...
                plan.s_hi-plan.s_lo);
            C_hat = evalMatrixChebBary(plan.validation_s(held),s_train, ...
                chebBarycentricWeights(q),C_train);
            error_by_q(isub,iq) = max(relativeActionError( ...
                C_hat,C_validation(:,:,held)));
        end

        qmax_local(isub) = error_by_q(isub,end);
        iq = find(error_by_q(isub,:) <= action_tolerance,1,'first');
        if isempty(iq)
            continue
        end

        q_local(isub) = q_candidates(iq);
        held_local(isub) = error_by_q(isub,iq);
        s_train = chebLobattoNodes(plan.s_lo,plan.s_hi,q_local(isub));
        C_train = lookupMatrices(to_delta(s_train),delta_pool,C_pool);
        C_random = lookupMatrices(to_delta(plan.random_s),delta_pool,C_pool);
        C_hat_random = evalMatrixChebBary(plan.random_s,s_train, ...
            chebBarycentricWeights(q_local(isub)),C_train);
        random_local(isub) = max(relativeActionError( ...
            C_hat_random,C_random));
    end

    retained = arrayfun(@(x) x.n_image(1),info_pool);
    results(in).Nclust = opt.Nclust;
    results(in).pass = all(isfinite(q_local)) && ...
        max(random_local) <= action_tolerance;
    results(in).q_local = q_local;
    results(in).q_max = max(q_local,[],'omitnan');
    results(in).held_out_error = max(held_local,[],'omitnan');
    results(in).random_error = max(random_local,[],'omitnan');
    results(in).qmax_test_error = max(qmax_local);
    results(in).retained_ellipse_nodes = [min(retained) max(retained)];
    results(in).ellipse_count_changes = sum(diff(retained)~=0);
    results(in).local_error_by_q = error_by_q;

    fprintf(['  q_local=%s; q=17 worst error %.3e; selected held-out ', ...
        '%.3e; random %.3e; retained nodes %d--%d; %s\n'], ...
        mat2str(q_local),results(in).qmax_test_error, ...
        results(in).held_out_error,results(in).random_error, ...
        results(in).retained_ellipse_nodes(1), ...
        results(in).retained_ellipse_nodes(2), ...
        passText(results(in).pass));
end

%% Summary and visible plots
fprintf('\n%-8s %-28s %6s %11s %11s %11s %14s\n', ...
    'Nclust','local q','q_max','q=17 error','held-out','random', ...
    'retained nodes');
for k = 1:numel(results)
    r = results(k);
    fprintf('%-8d %-28s %6s %11.3e %11s %11s %6d--%-6d\n', ...
        r.Nclust,mat2str(r.q_local),numberOrDash(r.q_max,'%d'), ...
        r.qmax_test_error,numberOrDash(r.held_out_error,'%.3e'), ...
        numberOrDash(r.random_error,'%.3e'), ...
        r.retained_ellipse_nodes(1),r.retained_ellipse_nodes(2));
end

f = figure('Name','Ellipse node density and local Cmap interpolation', ...
    'Color','w','Visible','on');
subplot(1,3,1);
plot([results.Nclust],[results.q_max],'-o','LineWidth',1.2);
ylabel('maximum selected local q'); xlabel('Nclust'); grid on;
title('Local interpolation length');

subplot(1,3,2);
semilogy([results.Nclust],[results.qmax_test_error],'-o', ...
    'LineWidth',1.2,'DisplayName','q=17 held-out'); hold on;
semilogy([results.Nclust],[results.random_error],'-s', ...
    'LineWidth',1.2,'DisplayName','selected-q random');
yline(action_tolerance,'--','DisplayName','target');
xlabel('Nclust'); ylabel('relative action error'); grid on;
legend('Location','best'); title('Accuracy');

subplot(1,3,3);
node_ranges = vertcat(results.retained_ellipse_nodes);
plot([results.Nclust],node_ranges(:,1),'-o','LineWidth',1.2, ...
    'DisplayName','minimum'); hold on;
plot([results.Nclust],node_ranges(:,2),'-s','LineWidth',1.2, ...
    'DisplayName','maximum');
xlabel('Nclust'); ylabel('retained ellipse sources per body'); grid on;
legend('Location','best'); title('Actual source counts');

saveas(f,figure_file);
drawnow;

save(results_file,'results','Nclust_values','n_subdomains', ...
    'q_candidates','q_validation','action_tolerance', ...
    'n_random_per_subdomain','rng_seed','plans','delta_star', ...
    'transition_margin','opt_base','R','-v7.3');
fprintf('\nResults saved to %s\n',relativePath(results_file,repo_root));
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
