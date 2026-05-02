function [chunk_pairs,opt] = resolveMobilityChunkPairs(total_pairs,opt,pool_size)
%RESOLVEMOBILITYCHUNKPAIRS Normalize the shared mobility chunk size.

if nargin < 1 || isempty(total_pairs)
    total_pairs = 0;
end
if nargin < 2 || isempty(opt)
    opt = struct();
end
if nargin < 3
    pool_size = [];
end

parallel_chunk = normalizeChunkValue( ...
    getOptField(opt,'parallel_precomp_chunk_pairs',[]), ...
    'parallel_precomp_chunk_pairs');
big_sparse_chunk = normalizeChunkValue( ...
    getOptField(opt,'mob_big_sparse_chunk_pairs',[]), ...
    'mob_big_sparse_chunk_pairs');

if ~isempty(parallel_chunk) && ~isempty(big_sparse_chunk) && ...
        parallel_chunk ~= big_sparse_chunk
    error('resolveMobilityChunkPairs:ConflictingChunkPairs', ...
        ['opt.parallel_precomp_chunk_pairs (%d) and ', ...
         'opt.mob_big_sparse_chunk_pairs (%d) must match.'], ...
        parallel_chunk,big_sparse_chunk);
end

if ~isempty(parallel_chunk)
    chunk_pairs = parallel_chunk;
elseif ~isempty(big_sparse_chunk)
    chunk_pairs = big_sparse_chunk;
else
    chunk_pairs = autoChunkPairs(total_pairs,pool_size);
end

if total_pairs > 0
    chunk_pairs = min(total_pairs,max(1,chunk_pairs));
else
    chunk_pairs = max(1,chunk_pairs);
end
chunk_pairs = double(chunk_pairs);

opt.parallel_precomp_chunk_pairs = chunk_pairs;
opt.mob_big_sparse_chunk_pairs = chunk_pairs;
end

function chunk_value = normalizeChunkValue(chunk_value,field_name)
if isempty(chunk_value)
    chunk_value = [];
    return
end

if ~(isnumeric(chunk_value) && isscalar(chunk_value) && ...
        isfinite(chunk_value) && chunk_value == round(chunk_value))
    error('resolveMobilityChunkPairs:BadChunkPairs', ...
        'opt.%s must be empty or an integer scalar.',field_name);
end

if chunk_value <= 0
    chunk_value = [];
else
    chunk_value = double(chunk_value);
end
end

function chunk_pairs = autoChunkPairs(total_pairs,pool_size)
pool_size = resolveAutoChunkPoolSize(pool_size);
chunk_pairs = max(1,ceil(max(1,total_pairs) / (4*pool_size)));
end

function pool_size = resolveAutoChunkPoolSize(pool_size)
if nargin < 1 || isempty(pool_size)
    pool_size = 0;
    if exist('gcp','file') == 2
        pool = gcp('nocreate');
        if ~isempty(pool)
            pool_size = pool.NumWorkers;
        end
    end
end

if ~(isnumeric(pool_size) && isscalar(pool_size) && isfinite(pool_size))
    pool_size = 0;
end
pool_size = max(1,round(pool_size));
end
