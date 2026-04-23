% Compare looped and batched pair_qv block updates in Laplace peanut matvecs.
%
% The benchmark isolates the projected-solve update
%
%   res(block) = res(block) + pair_qv_nonp(k)
%
% over all bodies touched by close-pair interactions. It compares the
% original per-body loop against the in-place batched
% reshape-and-column update used in the Laplace peanut matvecs.

script_name = mfilename;
script_date = 'Apr 23, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s (%s) ===\n',script_name,script_date);

if ~exist('P_sweep','var') || isempty(P_sweep)
    P_sweep = [10 40 100 250 1000];
end
if ~exist('R','var') || isempty(R)
    R = 2;
end
if ~exist('N_c','var') || isempty(N_c)
    N_c = 40;
end
if ~exist('N_f','var') || isempty(N_f)
    N_f = 80;
end
if ~exist('n_repeats','var') || isempty(n_repeats)
    n_repeats = 3;
end
if ~exist('rng_seed','var') || isempty(rng_seed)
    rng_seed = 230423;
end
if ~exist('pair_gap','var') || isempty(pair_gap)
    pair_gap = 0.02;
end

opt = getLaplace2Dparams(max(P_sweep),R,N_c,N_f);
a_c = opt.a_c;
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1).';
rbase_out_c = R*(cos(tout_c) + 1i*sin(tout_c));
N_large = numel(rbase_out_c);

rng(rng_seed,'twister');

fprintf('Pair-induced charge block update benchmark\n');
fprintf('gap=%g, N_c=%d, N_check=%d, repeats=%d\n', ...
    pair_gap,N_c,N_large,n_repeats);
fprintf('%8s %10s %10s %12s %12s %10s %12s\n', ...
    'P','n_pairs','n_touch','loop_s','batched_s','speedup','rel_inf');

for ip = 1:numel(P_sweep)
    P = P_sweep(ip);
    q = build_pair_test_centers(P,R,pair_gap);
    opt_pairs = opt;
    opt_pairs.P = P;
    [~,~,~,~,~,pairs] = getEnhancedGrid(q,opt_pairs);
    if isempty(pairs)
        error('apr23_laplace_pair_qv_block_update_compare:NoPairs', ...
            ['No close pairs found for P=%d with gap %.3g. ', ...
             'Reduce pair_gap to benchmark the pair_qv update.'], ...
            P,pair_gap);
    end

    has_neigh = false(1,P);
    has_neigh(pairs(:)) = true;

    pair_qv_nonp = randn(P,1);
    res0 = randn(P*N_large,1);

    loop_result = apply_pair_qv_loop(res0,pair_qv_nonp,pairs,N_large);
    batch_result = apply_pair_qv_batched(res0,pair_qv_nonp,pairs,P,N_large);
    rel_inf = norm(batch_result-loop_result,inf) / ...
        max(1,norm(loop_result,inf));

    t_loop = zeros(n_repeats,1);
    t_batch = zeros(n_repeats,1);
    for rep = 1:n_repeats
        t_loop(rep) = timeit(@() apply_pair_qv_loop( ...
            res0,pair_qv_nonp,pairs,N_large));
        t_batch(rep) = timeit(@() apply_pair_qv_batched( ...
            res0,pair_qv_nonp,pairs,P,N_large));
    end

    loop_s = median(t_loop);
    batch_s = median(t_batch);
    speedup = loop_s / max(batch_s,eps);
    fprintf('%8d %10d %10d %12.4e %12.4e %10.2f %12.3e\n', ...
        P,size(pairs,1),nnz(has_neigh),loop_s,batch_s,speedup,rel_inf);
end

function res = apply_pair_qv_loop(res,pair_qv_nonp,pairs,N_large)
has_neigh = sort(unique(pairs(:)));
for ii = 1:numel(has_neigh)
    k = has_neigh(ii);
    block = (k-1)*N_large+1:k*N_large;
    res(block) = res(block) + pair_qv_nonp(k);
end
end

function res = apply_pair_qv_batched(res,pair_qv_nonp,pairs,P,N_large)
has_neigh = false(1,P);
has_neigh(pairs(:)) = true;
if any(has_neigh)
    res_blocks = reshape(res,N_large,P);
    res_blocks(:,has_neigh) = bsxfun(@plus, ...
        res_blocks(:,has_neigh),reshape(pair_qv_nonp(has_neigh),1,[]));
    res = res_blocks(:);
end
end

function q = build_pair_test_centers(P,R,pair_gap)
q = (0:P-1).'*(2*R + pair_gap);
end
