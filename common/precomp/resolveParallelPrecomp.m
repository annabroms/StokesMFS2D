function [use_parallel,pool_size,pool_type] = resolveParallelPrecomp( ...
    requested_parallel,caller_name)
%RESOLVEPARALLELPRECOMP Resolve and, if needed, open the parallel pool.

if nargin < 1 || isempty(requested_parallel)
    requested_parallel = false;
end
if nargin < 2 || isempty(caller_name)
    caller_name = mfilename;
end

use_parallel = logical(requested_parallel);
pool_size = 0;
pool_type = 'none';
if ~use_parallel
    return
end

if isempty(ver('parallel')) || ...
        ~license('test','Distrib_Computing_Toolbox') || ...
        exist('gcp','file') ~= 2
    error([caller_name ':ParallelToolboxRequired'], ...
        ['opt.parallel_precomp or opt.parallel_big_sparse_build ', ...
         'requires Parallel Computing Toolbox. Open a pool before ', ...
         'benchmarking warm-pool speedups if you want to exclude ', ...
         'startup overhead.']);
end

pool = gcp('nocreate');
if isempty(pool)
    pool = gcp();
end

pool_size = pool.NumWorkers;
pool_type = classifyParallelPool(pool);
end

function pool_type = classifyParallelPool(pool)
pool_type = 'none';
if isempty(pool)
    return
end

pool_class = class(pool);
if strcmp(pool_class,'parallel.ProcessPool')
    pool_type = 'processes';
elseif strcmp(pool_class,'parallel.ThreadPool')
    pool_type = 'threads';
else
    pool_type = pool_class;
end
end
