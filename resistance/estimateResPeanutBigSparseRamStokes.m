function estimate = estimateResPeanutBigSparseRamStokes(P,N_c,N_check,n_pairs,opt)
%ESTIMATERESPEANUTBIGSPARSERAMSTOKES Estimate resistance sparse-map memory.
%
% opt.res_sparse_map_coarse selects the resistance coarse
% source-correction path:
%   false (default): keep the legacy pair_tau + source_scatter_rows
%                    update.
%   true:            build a direct global M_source_corr in addition to
%                    the pair-local maps that are still needed for u_corr
%                    and M_ft_corr.
% opt.res_smaller_ucorr replaces the combined sparse M_u_corr by a
% structured dense pair payload that reuses the shared one-body self block
% basis.Nii. The legacy aliases opt.sparse_map_coarse and
% opt.smaller_ucorr are still accepted.

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

chunk_pairs = max(1,round(getOptField(opt,'res_big_sparse_chunk_pairs',8)));
chunk_pairs = min(max(1,n_pairs),chunk_pairs);
u_corr_mode = resolveResBigSparseUCorrMode(opt);
sparse_map_coarse = logical(getOptField(opt,'res_sparse_map_coarse',...
    getOptField(opt,'sparse_map_coarse',false)));
combined_u_corr = strcmp(u_corr_mode,'combined');
smaller_u_corr = combined_u_corr && logical(getOptField(opt, ...
    'res_smaller_ucorr',getOptField(opt,'smaller_ucorr',false)));
combined_sparse_u_corr = combined_u_corr && ~smaller_u_corr;
if strcmp(u_corr_mode,'factored')
    velocity_correction = 'factored_sparse';
elseif smaller_u_corr
    velocity_correction = 'combined_structured';
else
    velocity_correction = 'combined_sparse';
end

n_coarse = P*N_c;
pair_rows = 4*N_c;
pair_total_rows = n_pairs*pair_rows;

local_pair = pair_rows*(4*N_c);
local_u = (4*N_check)*(4*N_c);
local_u_small = 2*(2*N_check)*(2*N_c);
local_ft = 6*(4*N_c);

counts = struct();
counts.source_corr = double(sparse_map_coarse)*n_pairs*local_pair;
counts.pair = n_pairs*local_pair;
counts.u = double(combined_sparse_u_corr)*n_pairs*2*local_u;
counts.u_small = double(smaller_u_corr)*n_pairs*local_u_small;
counts.u_cross = double(strcmp(u_corr_mode,'factored'))*n_pairs*local_u;
counts.u_peanut = double(strcmp(u_corr_mode,'factored'))*n_pairs*local_u;
counts.ft = n_pairs*local_ft;
counts.source_scatter_indices = double(~sparse_map_coarse)*n_pairs*(4*N_c);
counts.velocity_scatter_indices = n_pairs*(4*N_check);
counts.ft_scatter_indices = n_pairs*6;
counts.total = counts.source_corr + counts.pair + counts.u + ...
    counts.u_small + counts.u_cross + counts.u_peanut + counts.ft;

n_source_cols = 2*n_coarse;
sparse_cols = double(sparse_map_coarse)*n_source_cols + n_source_cols + ...
    double(combined_sparse_u_corr)*(n_source_cols + pair_total_rows) + ...
    double(strcmp(u_corr_mode,'factored'))*(n_source_cols + pair_total_rows) + ...
    n_source_cols;
chunk_entries = chunk_pairs*(local_pair + local_ft);
if sparse_map_coarse
    chunk_entries = chunk_entries + chunk_pairs*local_pair;
end
if combined_sparse_u_corr
    chunk_entries = chunk_entries + chunk_pairs*2*local_u;
elseif strcmp(u_corr_mode,'factored')
    chunk_entries = chunk_entries + chunk_pairs*2*local_u;
end

sparse_nnz = counts.source_corr + counts.pair + counts.u + ...
    counts.u_cross + counts.u_peanut + counts.ft;
dense_entries = counts.u_small;

estimate = struct();
estimate.matrix_plan = struct( ...
    'source_correction',ternary(sparse_map_coarse, ...
        'direct_sparse_global','direct_sparse_scatter'), ...
    'sparse_map_coarse',sparse_map_coarse, ...
    'velocity_correction',velocity_correction, ...
    'u_corr_mode',u_corr_mode, ...
    'combined_u_corr',combined_u_corr, ...
    'smaller_u_corr',smaller_u_corr, ...
    'force_torque_correction','direct_sparse', ...
    'chunk_pairs',chunk_pairs);
estimate.counts = counts;
estimate.big_sparse_matrix_bytes = 16*sparse_nnz + 8*sparse_cols + ...
    8*dense_entries;
estimate.big_sparse_auxiliary_bytes = 8*(counts.source_scatter_indices + ...
    counts.velocity_scatter_indices + counts.ft_scatter_indices);
estimate.big_sparse_build_bytes = 24*chunk_entries;
estimate.big_sparse_peak_bytes = estimate.big_sparse_matrix_bytes + ...
    estimate.big_sparse_auxiliary_bytes + estimate.big_sparse_build_bytes;
end

function out = ternary(cond,true_value,false_value)
if cond
    out = true_value;
else
    out = false_value;
end
end
