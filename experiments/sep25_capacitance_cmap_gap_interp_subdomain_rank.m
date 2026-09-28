%% SEP25_CAPACITANCE_CMAP_GAP_INTERP_SUBDOMAIN_RANK
% Test whether separate spatial bases on equal alpha subintervals reduce
% the rank of the gap-dependent Laplace capacitance pair map.
%
% The exact maps use Tikhonov regularisation, a constant ellipse
% discretisation with 150 candidate nodes, and no Cmap compression.  The
% two base gap ranges are each divided into 1, 2, 4, or 8 equal intervals
% in alpha=acosh(1+delta/(2R)).  A common left/right basis is built from
% q_train snapshots on every numerical interval and validated at the
% held-out nodes of a q_validation grid plus independent random gaps.
%
% reference_locations is deliberately exposed near the top.  "minimum"
% uses the smallest physical delta in each numerical interval, "midpoint"
% uses its arithmetic delta midpoint, and "maximum" is also supported.
% The experiment compares minimum and midpoint by default.
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

fprintf('=== sep25_capacitance_cmap_gap_interp_subdomain_rank ===\n\n');

%% Configuration: change reference_locations here
R = 2;
P = 2;
interpolation_coordinate = 'alpha';
subdomain_counts = [1 2 4 8];
reference_locations = {'minimum','midpoint'}; % also supports 'maximum'
primary_reference_location = 'midpoint';

q_train = 9;
q_validation = 33;
n_random_per_subdomain = 4;
spatial_action_tolerance = 1e-6;
snapshot_rank_tolerance = 1e-8;
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

switch interpolation_coordinate
    case 'alpha'
        coordinate_name = 'alpha = acosh(1+delta/(2R))';
        to_param = @(delta) acosh(1+delta/(2*R));
        to_delta = @(alpha) 2*R*(cosh(alpha)-1);
    case 'log_delta'
        coordinate_name = 's = log(delta/R)';
        to_param = @(delta) log(delta/R);
        to_delta = @(s) R*exp(s);
    otherwise
        error('Unknown interpolation coordinate "%s".', ...
            interpolation_coordinate);
end

base_ranges = struct( ...
    'name',{'small-gap range','large-gap range'}, ...
    'delta_lo',{delta_min,delta_star}, ...
    'delta_hi',{delta_star,delta_max});

assert(any(strcmp(reference_locations,primary_reference_location)), ...
    'primary_reference_location must occur in reference_locations.');
assert(mod(q_validation-1,q_train-1)==0, ...
    'The q_train grid must be nested in the q_validation grid.');

cache_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_subdomain_rank_cache.mat');
results_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_subdomain_rank_results.mat');
figure_file = fullfile(repo_root,'data', ...
    'sep25_capacitance_cmap_gap_interp_subdomain_rank.png');

fprintf(['Coordinate: %s; frozen ellipse Nclust=%d at delta_min/R=%.3g.\n' ...
    'q_train=%d, q_validation=%d, spatial action tolerance=%.1e.\n' ...
    'Reference locations: %s.\n\n'],coordinate_name,opt.Nclust, ...
    opt.smallest_delta/R,q_train,q_validation, ...
    spatial_action_tolerance,strjoin(reference_locations,', '));

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

%% Plan equal subintervals and all exact snapshots
rng(rng_seed);
plans = struct([]);
n_plans = numel(base_ranges)*sum(subdomain_counts);
n_delta_per_plan = q_train+q_validation+n_random_per_subdomain+ ...
    numel(reference_locations);
all_delta = zeros(n_plans*n_delta_per_plan,1);
delta_index = 0;
for ib = 1:numel(base_ranges)
    base = base_ranges(ib);
    x_base = [to_param(base.delta_lo),to_param(base.delta_hi)];
    for ic = 1:numel(subdomain_counts)
        n_sub = subdomain_counts(ic);
        edges = linspace(x_base(1),x_base(2),n_sub+1);
        for isub = 1:n_sub
            plan.base_index = ib;
            plan.n_subdomains = n_sub;
            plan.local_index = isub;
            plan.x_lo = edges(isub);
            plan.x_hi = edges(isub+1);
            plan.delta_lo = to_delta(plan.x_lo);
            plan.delta_hi = to_delta(plan.x_hi);
            plan.train_x = chebLobattoNodes( ...
                plan.x_lo,plan.x_hi,q_train);
            plan.validation_x = chebLobattoNodes( ...
                plan.x_lo,plan.x_hi,q_validation);
            plan.random_x = plan.x_lo+(plan.x_hi-plan.x_lo)* ...
                rand(n_random_per_subdomain,1);
            plans = [plans; plan]; %#ok<AGROW>

            new_delta = [to_delta(plan.train_x); ...
                to_delta(plan.validation_x); to_delta(plan.random_x)];
            indices = delta_index+(1:numel(new_delta));
            all_delta(indices) = new_delta;
            delta_index = delta_index+numel(new_delta);
            for irule = 1:numel(reference_locations)
                delta_index = delta_index+1;
                all_delta(delta_index) = referenceDelta( ...
                    plan.delta_lo,plan.delta_hi, ...
                    reference_locations{irule});
            end
        end
    end
end
all_delta = all_delta(1:delta_index);

delta_pool = mergeNearlyEqual(all_delta,1e-12);
fprintf('Building/reusing %d distinct exact Cmap snapshots...\n', ...
    numel(delta_pool));
[C_pool,~] = getOrBuildCmapSnapshots( ...
    delta_pool,R,opt,grids,cache_file);

%% Build and validate one spatial basis per numerical subinterval
record_template = struct('base_index',[],'base_name','','n_subdomains',[], ...
    'local_index',[],'delta_lo',[],'delta_hi',[], ...
    'reference_location','','reference_delta',[], ...
    'r_selected',[],'projection_action_error',[], ...
    'snapshot_rank',[],'rank_candidates',[],'candidate_errors',[]);
records = repmat(record_template,0,1);

for ip = 1:numel(plans)
    plan = plans(ip);
    C_train = lookupMatrices(to_delta(plan.train_x),delta_pool,C_pool);
    C_validation = lookupMatrices( ...
        to_delta(plan.validation_x),delta_pool,C_pool);
    held = heldOutMask(plan.validation_x,plan.train_x, ...
        plan.x_hi-plan.x_lo);
    C_random = lookupMatrices(to_delta(plan.random_x),delta_pool,C_pool);
    C_test = cat(3,C_validation(:,:,held),C_random);

    for irule = 1:numel(reference_locations)
        reference_location = reference_locations{irule};
        delta_ref = referenceDelta( ...
            plan.delta_lo,plan.delta_hi,reference_location);
        C_ref = lookupMatrices(delta_ref,delta_pool,C_pool);
        C_ref = C_ref(:,:,1);

        dC_train = C_train-C_ref;
        dC_test = C_test-C_ref;
        basis = commonBasisFromSnapshots(dC_train,r_max);
        [r_selected,selected_error,rank_candidates,candidate_errors] = ...
            selectSpatialRank(dC_test,C_test,basis.U,basis.V, ...
            spatial_action_tolerance,rank_search_stride);

        n = size(C_train,1);
        F = reshape(dC_train,n*n,q_train).';
        sigma_F = svd(F,'econ');
        snapshot_rank = sum(sigma_F > ...
            snapshot_rank_tolerance*max(sigma_F(1),eps));

        row = record_template;
        row.base_index = plan.base_index;
        row.base_name = base_ranges(plan.base_index).name;
        row.n_subdomains = plan.n_subdomains;
        row.local_index = plan.local_index;
        row.delta_lo = plan.delta_lo;
        row.delta_hi = plan.delta_hi;
        row.reference_location = reference_location;
        row.reference_delta = delta_ref;
        row.r_selected = r_selected;
        row.projection_action_error = selected_error;
        row.snapshot_rank = snapshot_rank;
        row.rank_candidates = rank_candidates;
        row.candidate_errors = candidate_errors;
        records(end+1,1) = row; %#ok<SAGROW>
    end
end

%% Aggregate and report full-range versus local ranks
summary_template = struct('base_index',[],'base_name','', ...
    'n_subdomains',[],'reference_location','','local_ranks',[], ...
    'max_rank',[],'max_projection_error',[],'snapshot_ranks',[]);
summary = repmat(summary_template,0,1);
for ib = 1:numel(base_ranges)
    for ic = 1:numel(subdomain_counts)
        n_sub = subdomain_counts(ic);
        for irule = 1:numel(reference_locations)
            rule = reference_locations{irule};
            take = find([records.base_index]==ib & ...
                [records.n_subdomains]==n_sub & ...
                strcmp({records.reference_location},rule));
            [~,order] = sort([records(take).local_index]);
            take = take(order);
            row = summary_template;
            row.base_index = ib;
            row.base_name = base_ranges(ib).name;
            row.n_subdomains = n_sub;
            row.reference_location = rule;
            row.local_ranks = [records(take).r_selected];
            row.max_rank = max(row.local_ranks);
            row.max_projection_error = max( ...
                [records(take).projection_action_error]);
            row.snapshot_ranks = [records(take).snapshot_rank];
            summary(end+1,1) = row; %#ok<SAGROW>
        end
    end
end

fprintf('\n=== Spatial rank by equal %s subintervals ===\n', ...
    interpolation_coordinate);
fprintf('%-17s %6s %-9s %-28s %7s %11s %-20s\n', ...
    'base range','pieces','reference','local spatial ranks','r_max', ...
    'test error','snapshot ranks');
for k = 1:numel(summary)
    row = summary(k);
    fprintf('%-17s %6d %-9s %-28s %7d %11.3e %-20s\n', ...
        row.base_name,row.n_subdomains,row.reference_location, ...
        mat2str(row.local_ranks),row.max_rank, ...
        row.max_projection_error,mat2str(row.snapshot_ranks));
end

fprintf('\n=== Primary reference: %s ===\n',primary_reference_location);
for ib = 1:numel(base_ranges)
    one = summary([summary.base_index]==ib & ...
        [summary.n_subdomains]==1 & ...
        strcmp({summary.reference_location},primary_reference_location));
    fprintf('%s: full-range r=%d.\n',base_ranges(ib).name,one.max_rank);
    for ic = 2:numel(subdomain_counts)
        row = summary([summary.base_index]==ib & ...
            [summary.n_subdomains]==subdomain_counts(ic) & ...
            strcmp({summary.reference_location}, ...
            primary_reference_location));
        fprintf('  %d pieces: local r=%s, maximum %d (change %+d).\n', ...
            row.n_subdomains,mat2str(row.local_ranks),row.max_rank, ...
            row.max_rank-one.max_rank);
    end
end

%% Explain the q=5 versus q=9 online arithmetic for the measured ranks
n = 2*opt.N_c;
primary_rows = summary(strcmp({summary.reference_location}, ...
    primary_reference_location));
measured_ranks = unique([primary_rows.local_ranks]);
fprintf('\n=== Online arithmetic: q=5 versus q=9 ===\n');
fprintf(['Full-matrix barycentric evaluation scales as q*n^2: q=9 uses ', ...
    '9/5=1.8 times the interpolation work of q=5.  Once an explicit ', ...
    'matrix is formed, its subsequent n-by-n matvec cost is independent ', ...
    'of q.\n']);
for r = measured_ranks
    cost5 = n^2+2*n*r+r^2+5*r^2;
    cost9 = n^2+2*n*r+r^2+9*r^2;
    fprintf(['r=%3d: C_ref plus reduced apply and core interpolation gives ', ...
        'cost(9)/cost(5)=%.3f in the leading-count model.\n'], ...
        r,cost9/cost5);
end

%% Visible figure
f = figure('Name','Per-subinterval Cmap spatial rank', ...
    'Color','w','Visible','on');
for ib = 1:numel(base_ranges)
    subplot(1,numel(base_ranges),ib);
    hold on;
    for irule = 1:numel(reference_locations)
        rule = reference_locations{irule};
        take = find([summary.base_index]==ib & ...
            strcmp({summary.reference_location},rule));
        plot([summary(take).n_subdomains],[summary(take).max_rank], ...
            '-o','LineWidth',1.2,'DisplayName',rule);
    end
    xlabel('equal alpha subintervals');
    ylabel('maximum local spatial rank');
    title(base_ranges(ib).name);
    legend('Location','best');
    grid on;
end
saveas(f,figure_file);
drawnow;

save(results_file,'records','summary','plans','base_ranges', ...
    'subdomain_counts','reference_locations','primary_reference_location', ...
    'q_train','q_validation','n_random_per_subdomain', ...
    'spatial_action_tolerance','snapshot_rank_tolerance','r_max', ...
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

function [r_selected,selected_error,ranks_tested,errors_tested] = ...
        selectSpatialRank(dC_test,C_test,U,V,tolerance,stride)
r_max = size(U,2);
coarse_ranks = unique([1,stride:stride:r_max,r_max]);
coarse_errors = nan(size(coarse_ranks));
first_pass = [];
for k = 1:numel(coarse_ranks)
    coarse_errors(k) = maxProjectionActionError( ...
        dC_test,C_test,U,V,coarse_ranks(k));
    if coarse_errors(k) <= tolerance
        first_pass = k;
        break
    end
end

if isempty(first_pass)
    r_selected = NaN;
    selected_error = NaN;
    ranks_tested = coarse_ranks;
    errors_tested = coarse_errors;
    return
end

if first_pass==1
    refine_ranks = coarse_ranks(1);
else
    refine_ranks = (coarse_ranks(first_pass-1)+1):coarse_ranks(first_pass);
end
refine_errors = nan(size(refine_ranks));
for k = 1:numel(refine_ranks)
    refine_errors(k) = maxProjectionActionError( ...
        dC_test,C_test,U,V,refine_ranks(k));
end
local_pass = find(refine_errors <= tolerance,1,'first');
r_selected = refine_ranks(local_pass);
selected_error = refine_errors(local_pass);
ranks_tested = [coarse_ranks(1:first_pass),refine_ranks];
errors_tested = [coarse_errors(1:first_pass),refine_errors];
end

function error = maxProjectionActionError(dC,C,U,V,r)
U_r = U(:,1:r);
V_r = V(:,1:r);
error = 0;
for k = 1:size(C,3)
    projected = U_r*(U_r'*dC(:,:,k)*V_r)*V_r';
    current = norm(dC(:,:,k)-projected,2)/max(norm(C(:,:,k),2),eps);
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

function p = relativePath(filename,root)
prefix = [root filesep];
if strncmp(filename,prefix,numel(prefix))
    p = filename(numel(prefix)+1:end);
else
    p = filename;
end
end
