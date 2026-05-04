function [chunk_pairs, opt] = resolveMobilityChunkPairs(total_pairs, opt, pool_size)
%RESOLVEMOBILITYCHUNKPAIRS Normalize the shared mobility chunk size.

    if nargin < 1 || isempty(total_pairs); total_pairs = 0;        end
    if nargin < 2 || isempty(opt);         opt         = struct();  end
    if nargin < 3;                         pool_size   = [];        end

    pc  = normalizeChunkValue(getOptField(opt, 'parallel_precomp_chunk_pairs', []),...
                              'parallel_precomp_chunk_pairs');
    bsc = normalizeChunkValue(getOptField(opt, 'mob_big_sparse_chunk_pairs',   []),...
                              'mob_big_sparse_chunk_pairs');

    if ~isempty(pc) && ~isempty(bsc) && pc ~= bsc
        error('resolveMobilityChunkPairs:ConflictingChunkPairs',...
            ['opt.parallel_precomp_chunk_pairs (%d) and '...
             'opt.mob_big_sparse_chunk_pairs (%d) must match.'], pc, bsc);
    end

    % First non-empty wins; fall back to auto
    if     ~isempty(pc);  chunk_pairs = pc;
    elseif ~isempty(bsc); chunk_pairs = bsc;
    else;                 chunk_pairs = autoChunkPairs(total_pairs, pool_size);
    end

    % Clamp to [1, total_pairs] (skip upper clamp if total unknown)
    chunk_pairs = max(1, chunk_pairs);
    if total_pairs > 0
        chunk_pairs = min(total_pairs, chunk_pairs);
    end
    chunk_pairs = double(chunk_pairs);

    opt.parallel_precomp_chunk_pairs = chunk_pairs;
    opt.mob_big_sparse_chunk_pairs   = chunk_pairs;
end

% -------------------------------------------------------------------------

function chunk = normalizeChunkValue(chunk, field_name)
    if isempty(chunk); return; end   % no preference expressed

    if ~(isnumeric(chunk) && isscalar(chunk) &&...
         isfinite(chunk)  && chunk == round(chunk))
        error('resolveMobilityChunkPairs:BadChunkPairs',...
            'opt.%s must be empty or an integer scalar.', field_name);
    end

    chunk = double(chunk);
    if chunk <= 0; chunk = []; end   % non-positive → treat as unset
end

% -------------------------------------------------------------------------

function chunk = autoChunkPairs(total_pairs, pool_size)
    pool_size = resolvePoolSize(pool_size);
    chunk     = max(1, ceil(total_pairs / (4 * pool_size)));
end

% -------------------------------------------------------------------------

function pool_size = resolvePoolSize(pool_size)
    % Auto-detect from active pool if not provided
    if isempty(pool_size) && exist('gcp', 'file') == 2
        pool = gcp('nocreate');
        if ~isempty(pool)
            pool_size = pool.NumWorkers;
        end
    end

    % Validate and clamp — invalid or missing → default to 1
    if isnumeric(pool_size) && isscalar(pool_size) && isfinite(pool_size)
        pool_size = max(1, round(pool_size));
    else
        pool_size = 1;
    end
end