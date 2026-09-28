function model = buildLaplaceCmapInterpolator(opt,rbase_in_c, ...
        rbase_in_f,rout_base_f)
%BUILDLAPLACECMAPINTERPOLATOR Train one adaptive production alpha model.
%
% This is an explicit, potentially expensive training operation.  Normal
% scripts should call prepareLaplaceCmapInterpolation, which loads a saved
% compatible model and asks before invoking this routine on a cache miss.

[opt,mode] = configureLaplaceCapacitanceInterpolation(opt);
if strcmp(mode,'none')
    error('buildLaplaceCmapInterpolator:InterpolationDisabled', ...
        ['Set opt.use_interpolation to ''reduced_noconst'', ', ...
         '''reduced'', or ''full''.']);
end

if nargin < 2 || isempty(rbase_in_c)
    [rbase_in_c,rbase_in_f,rout_base_f] = defaultGrids(opt);
end
rbase_in_c = rbase_in_c(:);
rbase_in_f = rbase_in_f(:);
rout_base_f = rout_base_f(:);
if numel(rbase_in_c) ~= opt.N_cmap || ...
        numel(rbase_in_f) ~= opt.N_f || ...
        numel(rout_base_f) ~= ceil(opt.a_f*opt.N_f)
    error('buildLaplaceCmapInterpolator:GridSizeMismatch', ...
        'The supplied canonical grids do not match the configured sizes.');
end

signature = buildLaplaceCmapInterpolationSignature(opt,true);
signature_id = laplaceInterpolationSignatureId(signature);
exact_signature = buildLaplaceCmapInterpolationSignature(opt,false);
exact_signature = rmfield(exact_signature,'mode');

q_candidates = signature.node_candidates;
q_basis = max(q_candidates);
q_validation = signature.validation_nodes;
panel_count_candidates = signature.panel_count_candidates;
audit_fraction = [0.137 0.367 0.643 0.883];
R = opt.rad;
to_alpha = @(delta) acosh(1+delta/(2*R));
to_delta = @(alpha) 2*R*(cosh(alpha)-1);
base_ranges = makeBaseRanges(opt,to_alpha);

[candidate_plans,all_delta] = planCandidates(base_ranges, ...
    panel_count_candidates,q_candidates,q_basis,q_validation, ...
    audit_fraction,to_delta);
all_delta = mergeNearlyEqual(all_delta,1e-12);
grids = struct('rbase_in_c',rbase_in_c, ...
    'rbase_in_f',rbase_in_f,'rout_base_f',rout_base_f);
fprintf(['Training %s Laplace Cmap model: C tolerance %.3e, ', ...
    'C_Q tolerance %.3e.\n'],mode,opt.interpolation_tol, ...
    opt.charge_interpolation_tol);
fprintf(['Searching panel counts %s and node counts %s on %d ', ...
    'base range(s).\n'],mat2str(panel_count_candidates), ...
    mat2str(q_candidates),numel(base_ranges));
[delta_pool,C_pool,QV_pool,n_new,exact_cache_file] = ...
    getOrBuildExactSnapshots(all_delta,opt,grids,exact_signature);
fprintf('Exact snapshot pool: %d total, %d newly built; cache %s.\n', ...
    numel(delta_pool),n_new,exact_cache_file);

n_base = numel(base_ranges);
n_count = numel(panel_count_candidates);
candidate_results = cell(n_base,n_count);
for ib = 1:n_base
    for ic = 1:n_count
        plans = candidate_plans{ib,ic};
        selections = repmat(emptySelection(),numel(plans),1);
        pass = true;
        for ip = 1:numel(plans)
            selections(ip) = selectPanel(plans(ip),mode,q_candidates, ...
                opt.interpolation_tol,opt.charge_interpolation_tol, ...
                delta_pool,C_pool,QV_pool);
            pass = pass && selections(ip).pass;
        end
        candidate_results{ib,ic} = struct('pass',pass, ...
            'plans',plans,'selections',selections);
        printCandidate(base_ranges(ib).name,panel_count_candidates(ic), ...
            selections,pass,mode);
    end
end

[chosen_counts,chosen] = chooseConfiguration(candidate_results, ...
    panel_count_candidates,mode);
if isempty(chosen_counts)
    error('buildLaplaceCmapInterpolator:NoPassingConfiguration', ...
        ['No configuration met C tolerance %.3e and C_Q tolerance ', ...
         '%.3e with at most %d panels per base range and %d nodes.'], ...
        opt.interpolation_tol,opt.charge_interpolation_tol, ...
        max(panel_count_candidates),max(q_candidates));
end

panel_template = emptyModelPanel();
n_panels = sum(chosen_counts);
model_panels = repmat(panel_template,n_panels,1);
index = 0;
for ib = 1:n_base
    item = chosen{ib};
    for ip = 1:numel(item.plans)
        index = index+1;
        model_panels(index) = materializePanel(item.plans(ip), ...
            item.selections(ip),mode,delta_pool,C_pool,QV_pool);
    end
end

model = struct();
model.kind = 'laplace_capacitance_cmap_alpha';
model.version = 3;
model.mode = mode;
model.coordinate = 'alpha';
model.R = R;
model.N_cmap = opt.N_cmap;
model.delta_min = opt.smallest_delta;
model.delta_max = opt.delta_pair;
model.delta_star = (R-opt.Rp_f)^2/opt.Rp_f;
model.action_tolerance = opt.interpolation_tol;
model.charge_action_tolerance = opt.charge_interpolation_tol;
model.signature = signature;
model.signature_id = signature_id;
model.panels = model_panels;
model.n_panels = n_panels;
model.panel_counts_by_base = chosen_counts;
model.n_exact_snapshots = numel(delta_pool);
model.n_new_exact_snapshots = n_new;
model.exact_cache_file = exact_cache_file;
model.q_by_panel = [model_panels.q];
model.q_charge_by_panel = [model_panels.q_charge];
model.spatial_ranks = [model_panels.spatial_rank];
model.rbase_in_c = rbase_in_c;
model.rbase_in_f = rbase_in_f;
model.rout_base_f = rout_base_f;
model.training = struct('q_candidates',q_candidates, ...
    'panel_count_candidates',panel_count_candidates, ...
    'q_basis',q_basis,'q_validation',q_validation, ...
    'audit_fraction',audit_fraction, ...
    'max_C_error',max([model_panels.certified_C_error]), ...
    'max_CQ_error',max([model_panels.certified_CQ_error]));
model.profile = struct('ellipse_constant',opt.ellipse_constant, ...
    'Nclust',opt.Nclust,'use_tikhonov',opt.use_tikhonov, ...
    'tikhonov_tol',opt.tikhonov_tol, ...
    'compress_cmap',opt.compress_cmap, ...
    'reference_location',referenceLocation(mode), ...
    'basis_nodes',q_basis);
model.model_file = '';

fprintf('Selected panel counts by base range: %s.\n', ...
    mat2str(chosen_counts));
reportLaplaceCmapInterpolator(model,'newly trained (unsaved)');
end

function ranges = makeBaseRanges(opt,to_alpha)
delta_min = opt.smallest_delta;
delta_max = opt.delta_pair;
delta_star = (opt.rad-opt.Rp_f)^2/opt.Rp_f;
if delta_star > delta_min && delta_star < delta_max
    edges = [delta_min delta_star delta_max];
else
    edges = [delta_min delta_max];
end
n = numel(edges)-1;
ranges = repmat(struct('name','','delta_lo',[], ...
    'delta_hi',[],'alpha_lo',[],'alpha_hi',[]),n,1);
for k = 1:n
    ranges(k).name = sprintf('base range %d',k);
    ranges(k).delta_lo = edges(k);
    ranges(k).delta_hi = edges(k+1);
    alpha = to_alpha(edges(k:k+1));
    ranges(k).alpha_lo = alpha(1);
    ranges(k).alpha_hi = alpha(2);
end
end

function [all_plans,all_delta] = planCandidates(ranges,counts, ...
        q_candidates,q_basis,q_validation,audit_fraction,to_delta)
all_plans = cell(numel(ranges),numel(counts));
delta_chunks = cell(0,1);
for ib = 1:numel(ranges)
    for ic = 1:numel(counts)
        count = counts(ic);
        edges = linspace(ranges(ib).alpha_lo,ranges(ib).alpha_hi, ...
            count+1);
        plans = repmat(emptyPlan(),count,1);
        for ip = 1:count
            plan = emptyPlan();
            plan.base_index = ib;
            plan.local_index = ip;
            plan.alpha_lo = edges(ip);
            plan.alpha_hi = edges(ip+1);
            plan.delta_lo = to_delta(plan.alpha_lo);
            plan.delta_hi = to_delta(plan.alpha_hi);
            plan.reference_delta = (plan.delta_lo+plan.delta_hi)/2;
            plan.train_alpha = cell(numel(q_candidates),1);
            for iq = 1:numel(q_candidates)
                plan.train_alpha{iq} = chebLobatto( ...
                    plan.alpha_lo,plan.alpha_hi,q_candidates(iq));
                delta_chunks{end+1,1} = ...
                    to_delta(plan.train_alpha{iq}); %#ok<AGROW>
            end
            plan.basis_alpha = chebLobatto( ...
                plan.alpha_lo,plan.alpha_hi,q_basis);
            plan.validation_alpha = chebLobatto( ...
                plan.alpha_lo,plan.alpha_hi,q_validation);
            plan.audit_alpha = plan.alpha_lo+(plan.alpha_hi-plan.alpha_lo)* ...
                audit_fraction(:);
            delta_chunks{end+1,1} = to_delta(plan.basis_alpha); %#ok<AGROW>
            delta_chunks{end+1,1} = to_delta(plan.validation_alpha); %#ok<AGROW>
            delta_chunks{end+1,1} = to_delta(plan.audit_alpha); %#ok<AGROW>
            delta_chunks{end+1,1} = plan.reference_delta; %#ok<AGROW>
            plans(ip) = plan;
        end
        all_plans{ib,ic} = plans;
    end
end
all_delta = vertcat(delta_chunks{:});
end

function selection = selectPanel(plan,mode,q_candidates,C_tol,CQ_tol, ...
        delta_pool,C_pool,QV_pool)
Cref = lookupStack(plan.reference_delta,delta_pool,C_pool);
Cref = Cref(:,:,1);
C_basis = lookupStack(alphaToDelta(plan.basis_alpha),delta_pool,C_pool);
C_validation = lookupStack(alphaToDelta(plan.validation_alpha), ...
    delta_pool,C_pool);
C_audit = lookupStack(alphaToDelta(plan.audit_alpha),delta_pool,C_pool);
QV_validation = lookupStack(alphaToDelta(plan.validation_alpha), ...
    delta_pool,QV_pool);
QV_audit = lookupStack(alphaToDelta(plan.audit_alpha),delta_pool,QV_pool);

selection = emptySelection();
if strcmp(mode,'reduced_noconst')
    selection.Cref = zeros(size(Cref));
else
    selection.Cref = Cref;
end
if isReducedMode(mode)
    basis = buildLaplaceCmapSpatialBasis( ...
        C_basis-selection.Cref,size(Cref,1));
    selection.U_full = basis.U;
    selection.V_full = basis.V;
    spatial_exact = cat(3,C_validation,C_audit);
    minimum_rank = selectProjectionRank(spatial_exact,selection.Cref, ...
        basis.U,basis.V,C_tol);
    if isfinite(minimum_rank)
        for rank = minimum_rank:size(Cref,1)
            for iq = 1:numel(q_candidates)
                q = q_candidates(iq);
                [C_error,C_validation_error,C_audit_error] = ...
                    reducedInterpolationError(plan.train_alpha{iq},q, ...
                    rank,selection.Cref,basis.U,basis.V,plan, ...
                    C_validation,C_audit, ...
                    delta_pool,C_pool);
                if C_error <= C_tol
                    selection.rank = rank;
                    selection.q = q;
                    selection.C_error = C_error;
                    selection.C_validation_error = C_validation_error;
                    selection.C_audit_error = C_audit_error;
                    break
                end
            end
            if isfinite(selection.rank)
                break
            end
        end
    end
else
    for iq = 1:numel(q_candidates)
        q = q_candidates(iq);
        [C_error,C_validation_error,C_audit_error] = ...
            fullInterpolationError(plan.train_alpha{iq},q,plan, ...
            C_validation,C_audit,delta_pool,C_pool);
        if C_error <= C_tol
            selection.q = q;
            selection.C_error = C_error;
            selection.C_validation_error = C_validation_error;
            selection.C_audit_error = C_audit_error;
            break
        end
    end
end

for iq = 1:numel(q_candidates)
    q = q_candidates(iq);
    alpha_train = plan.train_alpha{iq};
    QV_train = lookupStack(alphaToDelta(alpha_train),delta_pool,QV_pool);
    held = heldOutMask(plan.validation_alpha,alpha_train, ...
        plan.alpha_hi-plan.alpha_lo);
    QV_hat_validation = evalBarycentric(plan.validation_alpha(held), ...
        alpha_train,barycentricWeights(q),QV_train);
    QV_hat_audit = evalBarycentric(plan.audit_alpha,alpha_train, ...
        barycentricWeights(q),QV_train);
    validation_error = maxRelativeError( ...
        QV_hat_validation,QV_validation(:,:,held));
    audit_error = maxRelativeError(QV_hat_audit,QV_audit);
    if max(validation_error,audit_error) <= CQ_tol
        selection.q_charge = q;
        selection.CQ_error = max(validation_error,audit_error);
        selection.CQ_validation_error = validation_error;
        selection.CQ_audit_error = audit_error;
        break
    end
end
selection.pass = isfinite(selection.q) && isfinite(selection.q_charge) && ...
    (~isReducedMode(mode) || isfinite(selection.rank));

    function delta = alphaToDelta(alpha)
        R_local = inferRadius(plan.alpha_lo,plan.delta_lo);
        delta = 2*R_local*(cosh(alpha)-1);
    end
end

function R = inferRadius(alpha,delta)
R = delta/(2*(cosh(alpha)-1));
end

function rank = selectProjectionRank(C_exact,Cref,U,V,tolerance)
n = size(Cref,1);
coarse = unique([1 4:4:n n]);
rank = NaN;
previous = 0;
for k = 1:numel(coarse)
    r = coarse(k);
    if projectionError(C_exact,Cref,U,V,r) <= tolerance
        for candidate = previous+1:r
            if projectionError(C_exact,Cref,U,V,candidate) <= tolerance
                rank = candidate;
                return
            end
        end
    end
    previous = r;
end
end

function value = projectionError(C_exact,Cref,U,V,rank)
U = U(:,1:rank);
V = V(:,1:rank);
value = 0;
for k = 1:size(C_exact,3)
    dC = C_exact(:,:,k)-Cref;
    approximation = Cref+U*(U'*dC*V)*V';
    value = max(value,relativeError(approximation,C_exact(:,:,k)));
end
end

function [total_error,validation_error,audit_error] = ...
        reducedInterpolationError(alpha_train,q,rank,Cref,U,V,plan, ...
        C_validation,C_audit,delta_pool,C_pool)
C_train = lookupStack(alphaToDeltaLocal(alpha_train,plan),delta_pool,C_pool);
held = heldOutMask(plan.validation_alpha,alpha_train, ...
    plan.alpha_hi-plan.alpha_lo);
C_hat_validation = evalBarycentric(plan.validation_alpha(held), ...
    alpha_train,barycentricWeights(q),C_train);
C_hat_audit = evalBarycentric(plan.audit_alpha,alpha_train, ...
    barycentricWeights(q),C_train);
U = U(:,1:rank);
V = V(:,1:rank);
C_hat_validation = reduceStack(C_hat_validation,Cref,U,V);
C_hat_audit = reduceStack(C_hat_audit,Cref,U,V);
validation_error = maxRelativeError( ...
    C_hat_validation,C_validation(:,:,held));
audit_error = maxRelativeError(C_hat_audit,C_audit);
total_error = max(validation_error,audit_error);
end

function [total_error,validation_error,audit_error] = ...
        fullInterpolationError(alpha_train,q,plan,C_validation,C_audit, ...
        delta_pool,C_pool)
C_train = lookupStack(alphaToDeltaLocal(alpha_train,plan),delta_pool,C_pool);
held = heldOutMask(plan.validation_alpha,alpha_train, ...
    plan.alpha_hi-plan.alpha_lo);
C_hat_validation = evalBarycentric(plan.validation_alpha(held), ...
    alpha_train,barycentricWeights(q),C_train);
C_hat_audit = evalBarycentric(plan.audit_alpha,alpha_train, ...
    barycentricWeights(q),C_train);
validation_error = maxRelativeError( ...
    C_hat_validation,C_validation(:,:,held));
audit_error = maxRelativeError(C_hat_audit,C_audit);
total_error = max(validation_error,audit_error);
end

function delta = alphaToDeltaLocal(alpha,plan)
R = inferRadius(plan.alpha_lo,plan.delta_lo);
delta = 2*R*(cosh(alpha)-1);
end

function C = reduceStack(C,Cref,U,V)
for k = 1:size(C,3)
    dC = C(:,:,k)-Cref;
    C(:,:,k) = Cref+U*(U'*dC*V)*V';
end
end

function [chosen_counts,chosen] = chooseConfiguration(results,counts,mode)
n_base = size(results,1);
n_count = numel(counts);
if n_base == 1
    combinations = (1:n_count)';
else
    [first,second] = ndgrid(1:n_count,1:n_count);
    combinations = [first(:) second(:)];
end
best_score = [];
best_indices = [];
for k = 1:size(combinations,1)
    selections = repmat(emptySelection(),0,1);
    ok = true;
    for ib = 1:n_base
        item = results{ib,combinations(k,ib)};
        ok = ok && item.pass;
        selections = [selections;item.selections]; %#ok<AGROW>
    end
    if ~ok
        continue
    end
    q = [selections.q];
    q_charge = [selections.q_charge];
    if isReducedMode(mode)
        ranks = [selections.rank];
        score = [max(ranks) sum(ranks) max(q) max(q_charge) ...
            sum(q)+sum(q_charge) numel(q)];
    else
        score = [max(q) max(q_charge) sum(q)+sum(q_charge) numel(q)];
    end
    if isempty(best_score) || lexicographicLess(score,best_score)
        best_score = score;
        best_indices = combinations(k,:);
    end
end
if isempty(best_indices)
    chosen_counts = [];
    chosen = {};
    return
end
chosen_counts = counts(best_indices);
chosen = cell(n_base,1);
for ib = 1:n_base
    chosen{ib} = results{ib,best_indices(ib)};
end
end

function yes = lexicographicLess(a,b)
index = find(a~=b,1,'first');
yes = ~isempty(index) && a(index)<b(index);
end

function panel = materializePanel(plan,selection,mode,delta_pool,C_pool,QV_pool)
panel = emptyModelPanel();
panel.base_index = plan.base_index;
panel.local_index = plan.local_index;
panel.alpha_lo = plan.alpha_lo;
panel.alpha_hi = plan.alpha_hi;
panel.delta_lo = plan.delta_lo;
panel.delta_hi = plan.delta_hi;
panel.q = selection.q;
panel.alpha_nodes = plan.train_alpha{findNodeIndex(plan,selection.q)};
panel.barycentric_weights = barycentricWeights(selection.q);
panel.q_charge = selection.q_charge;
panel.charge_alpha_nodes = ...
    plan.train_alpha{findNodeIndex(plan,selection.q_charge)};
panel.charge_barycentric_weights = ...
    barycentricWeights(selection.q_charge);
panel.QV_nodes = lookupStack(alphaToDeltaLocal( ...
    panel.charge_alpha_nodes,plan),delta_pool,QV_pool);
panel.reference_delta = plan.reference_delta;
panel.Cref = selection.Cref;
panel.certified_C_error = selection.C_error;
panel.certified_CQ_error = selection.CQ_error;
panel.validation_C_error = selection.C_validation_error;
panel.audit_C_error = selection.C_audit_error;
panel.validation_CQ_error = selection.CQ_validation_error;
panel.audit_CQ_error = selection.CQ_audit_error;
C_nodes = lookupStack(alphaToDeltaLocal(panel.alpha_nodes,plan), ...
    delta_pool,C_pool);
if strcmp(mode,'full')
    panel.C_nodes = C_nodes;
    panel.Cref = [];
else
    projection_reference = panel.Cref;
    panel.spatial_rank = selection.rank;
    panel.U = selection.U_full(:,1:selection.rank);
    panel.V = selection.V_full(:,1:selection.rank);
    panel.B_nodes = projectStack( ...
        C_nodes-projection_reference,panel.U,panel.V);
    panel.training_projection_error = projectionError( ...
        lookupStack(alphaToDeltaLocal(plan.validation_alpha,plan), ...
        delta_pool,C_pool),projection_reference,selection.U_full, ...
        selection.V_full,selection.rank);
    if strcmp(mode,'reduced_noconst')
        panel.Cref = [];
    end
end
end

function index = findNodeIndex(plan,q)
index = find(cellfun(@numel,plan.train_alpha)==q,1);
end

function printCandidate(name,count,selections,pass,mode)
if pass
    if isReducedMode(mode)
        rank_text = mat2str([selections.rank]);
    else
        rank_text = 'n/a';
    end
    fprintf(['  %-12s panels=%d: PASS ranks=%s, C nodes=%s, ', ...
        'C_Q nodes=%s\n'],name,count,rank_text, ...
        mat2str([selections.q]),mat2str([selections.q_charge]));
else
    fprintf('  %-12s panels=%d: no passing configuration\n',name,count);
end
end

function tf = isReducedMode(mode)
tf = any(strcmp(mode,{'reduced','reduced_noconst'}));
end

function location = referenceLocation(mode)
if strcmp(mode,'reduced_noconst')
    location = 'none';
else
    location = 'delta_midpoint';
end
end

function [delta_pool,C_pool,QV_pool,n_new,cache_file] = ...
        getOrBuildExactSnapshots(delta_query,opt,grids,signature)
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cache_dir = fullfile(repo_root,'data','laplace_cmap_interpolation');
signature_id = laplaceInterpolationSignatureId(signature);
cache_file = fullfile(cache_dir,['exact_' signature_id '.mat']);
delta_pool = zeros(0,1);
C_pool = zeros(2*opt.N_cmap,2*opt.N_cmap,0);
QV_pool = zeros(2,2*opt.N_cmap,0);
if isfile(cache_file)
    loaded = load(cache_file,'cache');
    if isfield(loaded,'cache') && ...
            isequaln(loaded.cache.signature,signature)
        delta_pool = loaded.cache.delta;
        C_pool = loaded.cache.C;
        QV_pool = loaded.cache.QV;
    end
end

n_new = 0;
for k = 1:numel(delta_query)
    delta = delta_query(k);
    if isempty(findDelta(delta,delta_pool))
        snapshot = buildCanonicalLaplacePairMap(delta,opt,grids);
        delta_pool(end+1,1) = delta; %#ok<AGROW>
        C_pool(:,:,end+1) = snapshot.Cmap; %#ok<AGROW>
        QV_pool(:,:,end+1) = snapshot.Cmap_QV; %#ok<AGROW>
        n_new = n_new+1;
        if logical(getOptField(opt,'show_counter',false)) && ...
                (mod(n_new,10)==0 || k==numel(delta_query))
            fprintf('  exact interpolation snapshot %d/%d new\n', ...
                n_new,numel(delta_query));
        end
    end
end
[delta_pool,order] = sort(delta_pool);
C_pool = C_pool(:,:,order);
QV_pool = QV_pool(:,:,order);
if n_new > 0
    if ~isfolder(cache_dir)
        mkdir(cache_dir);
    end
    cache = struct('signature',signature,'delta',delta_pool, ...
        'C',C_pool,'QV',QV_pool);
    save(cache_file,'cache','-v7.3');
end
end

function index = findDelta(delta,pool)
index = find(abs(pool-delta) <= 1e-10*max(1,abs(delta)),1);
end

function stack = lookupStack(query,pool,stack_pool)
query = query(:);
stack = zeros(size(stack_pool,1),size(stack_pool,2),numel(query));
for k = 1:numel(query)
    index = findDelta(query(k),pool);
    if isempty(index)
        error('buildLaplaceCmapInterpolator:MissingSnapshot', ...
            'A planned exact snapshot is absent from the pool.');
    end
    stack(:,:,k) = stack_pool(:,:,index);
end
end

function values = evalBarycentric(query,nodes,weights,snapshots)
values = zeros(size(snapshots,1),size(snapshots,2),numel(query));
for k = 1:numel(query)
    difference = query(k)-nodes;
    [distance,index] = min(abs(difference));
    if distance <= 1e-12*max(1,max(abs(nodes)))
        values(:,:,k) = snapshots(:,:,index);
    else
        terms = weights./difference;
        coefficients = terms/sum(terms);
        values(:,:,k) = reshape(reshape(snapshots,[],numel(nodes))* ...
            coefficients,size(snapshots,1),size(snapshots,2));
    end
end
end

function error_value = maxRelativeError(approximation,exact)
error_value = 0;
for k = 1:size(exact,3)
    error_value = max(error_value,relativeError( ...
        approximation(:,:,k),exact(:,:,k)));
end
end

function value = relativeError(approximation,exact)
value = norm(approximation-exact,2)/max(norm(exact,2),eps);
end

function mask = heldOutMask(validation_nodes,training_nodes,scale)
distance = min(abs(validation_nodes(:)-training_nodes(:)'),[],2);
mask = distance > 1e-10*max(1,abs(scale));
end

function B = projectStack(dC,U,V)
B = zeros(size(U,2),size(V,2),size(dC,3));
for k = 1:size(dC,3)
    B(:,:,k) = U'*dC(:,:,k)*V;
end
end

function [rbase_in_c,rbase_in_f,rout_base_f] = defaultGrids(opt)
tc = (0:opt.N_cmap-1)'*(2*pi/opt.N_cmap);
rbase_in_c = opt.Rp_c*exp(1i*tc);
tf = (0:opt.N_f-1)'*(2*pi/opt.N_f);
rbase_in_f = opt.Rp_f*exp(1i*tf);
nout_f = ceil(opt.a_f*opt.N_f);
to = (0:nout_f-1)'*(2*pi/nout_f);
rout_base_f = opt.rad*exp(1i*to);
end

function nodes = chebLobatto(lo,hi,q)
j = (0:q-1)';
nodes = (lo+hi)/2+(hi-lo)/2*cos(j*pi/(q-1));
end

function weights = barycentricWeights(q)
j = (0:q-1)';
weights = (-1).^j;
weights([1 end]) = weights([1 end])/2;
end

function values = mergeNearlyEqual(values,relative_tolerance)
values = sort(values(:));
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

function plan = emptyPlan()
plan = struct('base_index',[],'local_index',[],'alpha_lo',[], ...
    'alpha_hi',[],'delta_lo',[],'delta_hi',[], ...
    'reference_delta',[],'train_alpha',{{}},'basis_alpha',[], ...
    'validation_alpha',[],'audit_alpha',[]);
end

function selection = emptySelection()
selection = struct('pass',false,'rank',NaN,'q',NaN,'q_charge',NaN, ...
    'C_error',Inf,'CQ_error',Inf,'C_validation_error',Inf, ...
    'C_audit_error',Inf,'CQ_validation_error',Inf, ...
    'CQ_audit_error',Inf,'Cref',[],'U_full',[],'V_full',[]);
end

function panel = emptyModelPanel()
panel = struct('base_index',[],'local_index',[], ...
    'alpha_lo',[],'alpha_hi',[],'delta_lo',[],'delta_hi',[], ...
    'q',[],'alpha_nodes',[],'barycentric_weights',[], ...
    'q_charge',[],'charge_alpha_nodes',[], ...
    'charge_barycentric_weights',[],'spatial_rank',[], ...
    'C_nodes',[],'QV_nodes',[],'reference_delta',[],'Cref',[], ...
    'U',[],'V',[],'B_nodes',[],'training_projection_error',[], ...
    'certified_C_error',[],'certified_CQ_error',[], ...
    'validation_C_error',[],'audit_C_error',[], ...
    'validation_CQ_error',[],'audit_CQ_error',[]);
end
