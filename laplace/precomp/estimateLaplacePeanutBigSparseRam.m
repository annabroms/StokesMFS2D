function estimate = estimateLaplacePeanutBigSparseRam(P,N_c,N_check,n_pairs,opt)
%ESTIMATELAPLACEPEANUTBIGSPARSERAM Estimate Laplace big-sparse memory.

if nargin < 5 || isempty(opt)
    opt = struct();
end

validateattributes(P,{'numeric'},{'scalar','integer','nonnegative'}, ...
    mfilename,'P',1);
validateattributes(N_c,{'numeric'},{'scalar','integer','positive'}, ...
    mfilename,'N_c',2);
validateattributes(N_check,{'numeric'},{'scalar','integer','positive'}, ...
    mfilename,'N_check',3);
validateattributes(n_pairs,{'numeric'},{'scalar','integer','nonnegative'}, ...
    mfilename,'n_pairs',4);

build_mode = lower(char(getOptField(opt,'lap_big_sparse_build_mode','auto')));
switch build_mode
    case 'auto'
        if logical(getOptField(opt,'get_bndry_field',true))
            build_mode = 'precomputed';
        else
            build_mode = 'streaming';
        end
    case {'precomputed','streaming'}
        % supported
    otherwise
        error('estimateLaplacePeanutBigSparseRam:BadMode', ...
            ['opt.lap_big_sparse_build_mode must be ''auto'', ', ...
             '''precomputed'', or ''streaming''.']);
end
chunk_pairs = max(1,round(getOptField(opt,'lap_big_sparse_chunk_pairs',8)));
chunk_pairs = min(max(1,n_pairs),chunk_pairs);
use_factored_pair_maps = logical(getOptField(opt, ...
    'reuse_pair_basis_by_sep',false));

n_coarse = P*N_c;
pair_rows = 2*N_c;
pair_total_rows = n_pairs*pair_rows;
n_u_rows = P*N_check;
n_qv_rows = 2*n_pairs;

local_pair_nonp = pair_rows*(2*N_c);
local_u_cross = (2*N_check)*(2*N_c);
local_u_peanut = (2*N_check)*pair_rows;
local_u_qv = (2*N_check)*2;
local_qv = 2*(2*N_c);

counts = struct();
if use_factored_pair_maps
    counts.pair_nonp = 0;
    counts.pair_nonp_canon = n_pairs*local_pair_nonp;
    counts.qv = 0;
    counts.qv_canon = n_pairs*local_qv;
else
    counts.pair_nonp = n_pairs*local_pair_nonp;
    counts.pair_nonp_canon = 0;
    counts.qv = n_pairs*local_qv;
    counts.qv_canon = 0;
end
counts.u_cross = n_pairs*local_u_cross;
counts.u_peanut = n_pairs*local_u_peanut;
counts.u_qv = n_pairs*local_u_qv;
counts.u = counts.u_cross + counts.u_peanut + counts.u_qv;
counts.source_scatter_indices = n_pairs*(4*N_c);
counts.qv_scatter_indices = n_pairs*2;
counts.dense_pair_projector = pair_rows*pair_rows;
counts.total = counts.pair_nonp + counts.pair_nonp_canon + counts.u + ...
    counts.qv + counts.qv_canon;

local_pair_entries = local_pair_nonp + local_qv;
chunk_entries = chunk_pairs*(local_pair_entries + local_u_cross + ...
    local_u_peanut + local_u_qv);
if use_factored_pair_maps
    sparse_cols = pair_total_rows + n_coarse + pair_total_rows + ...
        max(1,n_qv_rows) + pair_total_rows;
else
    sparse_cols = n_coarse + n_coarse + pair_total_rows + ...
        max(1,n_qv_rows) + n_coarse;
end
auxiliary_bytes = 8*(counts.source_scatter_indices + ...
    counts.qv_scatter_indices + counts.dense_pair_projector);

estimate = struct();
estimate.matrix_plan = struct('source_correction','factored_structured', ...
    'field_correction','factored_sparse', ...
    'pair_map_mode',ternary(use_factored_pair_maps,'rotated_factored', ...
        'direct'), ...
    'build_mode',build_mode, ...
    'chunk_pairs',chunk_pairs);
estimate.counts = counts;
if use_factored_pair_maps
    estimate.sparse_shape = struct( ...
        'M_pair_nonp_canon',[pair_total_rows pair_total_rows], ...
        'M_u_cross',[n_u_rows n_coarse], ...
        'M_u_peanut',[n_u_rows pair_total_rows], ...
        'M_u_qv',[n_u_rows n_qv_rows], ...
        'M_pair_qv_canon',[n_qv_rows pair_total_rows]);
else
    estimate.sparse_shape = struct( ...
        'M_pair_nonp',[pair_total_rows n_coarse], ...
        'M_u_cross',[n_u_rows n_coarse], ...
        'M_u_peanut',[n_u_rows pair_total_rows], ...
        'M_u_qv',[n_u_rows n_qv_rows], ...
        'M_pair_qv',[n_qv_rows n_coarse]);
end
estimate.estimated_sparse_bytes = 16*counts.total + 8*sparse_cols;
estimate.estimated_auxiliary_bytes = auxiliary_bytes;
estimate.estimated_build_bytes = 24*chunk_entries;
estimate.estimated_peak_bytes = estimate.estimated_sparse_bytes + ...
    estimate.estimated_build_bytes + estimate.estimated_auxiliary_bytes;
estimate.estimated_sparse_MB = estimate.estimated_sparse_bytes/1024^2;
estimate.estimated_auxiliary_MB = estimate.estimated_auxiliary_bytes/1024^2;
estimate.estimated_build_MB = estimate.estimated_build_bytes/1024^2;
estimate.estimated_peak_MB = estimate.estimated_peak_bytes/1024^2;
end

function out = ternary(cond,a,b)
if cond
    out = a;
else
    out = b;
end
end
