function [C_all,info_all] = getOrBuildCmapSnapshots(delta_vals,R,opt,grids,cache_file)
%GETORBUILDCMAPSNAPSHOTS Exact C(delta) snapshots, cached across script runs.
%
% [C_all,info_all] = getOrBuildCmapSnapshots(delta_vals,R,opt,grids,cache_file)
% returns C_all (n x n x numel(delta_vals)) and info_all (struct array),
% one exact build per requested delta (see buildCanonicalPairCmap). Builds
% already present in cache_file (matched by delta, R, and the relevant opt
% fields, to a relative tolerance of 1e-10) are reused; only genuinely new
% delta values are computed, so refining a nested Chebyshev grid from
% q=9->17->33->65 only pays for the newly introduced points. cache_file is
% expected to live under data/ (git-ignored); this function creates its
% parent folder if needed.
%
% Anna Broms, Sep 2026

delta_vals = delta_vals(:);
n_req = numel(delta_vals);

tagged_opt = struct('R',R,'N_c',opt.N_c,'N_f',opt.N_f,'N_peanut',opt.N_peanut, ...
    'a_c',opt.a_c,'a_f',opt.a_f,'delta_pair',opt.delta_pair, ...
    'Rp_c',opt.Rp_c,'Rp_f',opt.Rp_f, ...
    'ellipse_constant',opt.ellipse_constant, ...
    'smallest_delta',getOptField(opt,'smallest_delta',[]), ...
    'beta',opt.beta,'Nclust',opt.Nclust, ...
    'use_tikhonov',logical(getOptField(opt,'use_tikhonov',false)), ...
    'tikhonov_tol',getOptField(opt,'tikhonov_tol',[]), ...
    'compress_cmap',logical(getOptField(opt,'compress_cmap',false)));

cache = struct('delta',{},'C',{},'info',{},'tag',{});
if exist(cache_file,'file')
    loaded = load(cache_file,'cache');
    if isfield(loaded,'cache')
        cache = loaded.cache;
    end
end

C_all = zeros(2*opt.N_c,2*opt.N_c,n_req);
info_all = repmat(struct('n',[],'rank_fine',[],'rank_peanut',[], ...
    'effective_dof_fine',[],'effective_dof_peanut',[], ...
    'n_image',[],'t_build',[]),n_req,1);
n_new = 0;

for k = 1:n_req
    delta = delta_vals(k);
    hit = find_cached(cache,delta,tagged_opt);
    if isempty(hit)
        [C,info] = buildCanonicalPairCmap(delta,R,opt,grids);
        entry.delta = delta;
        entry.C = C;
        entry.info = info;
        entry.tag = tagged_opt;
        cache(end+1) = entry; %#ok<AGROW>
        n_new = n_new+1;
    else
        C = cache(hit).C;
        info = cache(hit).info;
    end
    C_all(:,:,k) = C;
    info_all(k) = info;
end

if n_new > 0
    cache_dir = fileparts(cache_file);
    if ~isempty(cache_dir) && ~isfolder(cache_dir)
        mkdir(cache_dir);
    end
    save(cache_file,'cache','-v7.3');
end
fprintf('  getOrBuildCmapSnapshots: %d requested, %d newly built, %d reused from cache.\n', ...
    n_req,n_new,n_req-n_new);

end

function hit = find_cached(cache,delta,tag)
hit = [];
if isempty(cache)
    return
end
for i = 1:numel(cache)
    if ~isequal(cache(i).tag,tag)
        continue
    end
    if abs(cache(i).delta-delta) < 1e-10*max(1,abs(delta))
        hit = i;
        return
    end
end
end
