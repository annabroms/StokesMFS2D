function estimate = estimateMobPeanutBigSparseRamStokes(P,N_c,N_check,n_pairs,opt)
%ESTIMATEMOBPEANUTBIGSPARSERAMSTOKES Estimate big-sparse matrix memory.
%
% estimate = estimateMobPeanutBigSparseRamStokes(P,N_c,N_check,n_pairs,opt)
% returns conservative byte counts for the sparse matrices and temporary
% sparse-entry arrays used by buildMobPeanutBigSparseStokes. Auxiliary
% arrays cover dense local projectors and indexed source scatter rows used
% by opt.big_sparse_structured_apply.
%
% The estimate needs the number of close pairs. For a solve-grid setup this
% can be obtained with getEnhancedGrid or count_close_pairs before building
% the matrices.

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

plan = getMobPeanutBigSparseMatrixPlan(opt);

n_coarse = P*N_c;
n_source_cols = 2*n_coarse;
n_u_cols = 2*n_coarse;
pair_rows = 4*N_c;
pair_total_rows = n_pairs*pair_rows;

source_entries = n_pairs*(8*N_c)*(4*N_c);
u_entries = n_pairs*(4*N_check)*(4*N_c);
pair_nonp_entries = n_pairs*pair_rows*(4*N_c);
pair_projector_entries = n_pairs*pair_rows*pair_rows;
source_scatter_entries = n_pairs*(8*N_c);
source_scatter_index_entries = n_pairs*(8*N_c);
dense_pair_projector_entries = pair_rows*pair_rows;

counts = struct();
counts.source = double(plan.direct_source_corr)*source_entries;
counts.u = double(plan.need_u_corr_matrix)*u_entries;
counts.pair_nonp = double(plan.need_pair_nonp_matrix)*pair_nonp_entries;
counts.pair_nonp_u_corr = double(plan.combined_pair_u_corr)* ...
    (pair_nonp_entries + u_entries);
counts.pair_projector = double(plan.need_pair_projector)*pair_projector_entries;
counts.source_scatter = double(plan.need_source_scatter)*source_scatter_entries;
counts.source_scatter_indices = double(plan.need_source_scatter_indices)* ...
    source_scatter_index_entries;
counts.dense_pair_projector = double(plan.need_dense_pair_projector)* ...
    dense_pair_projector_entries;
counts.u_cross = double(plan.need_u_cross)*u_entries;
counts.u_peanut = double(plan.need_u_peanut)*u_entries;
counts.total = counts.source + counts.u + counts.pair_nonp + ...
    counts.pair_nonp_u_corr + counts.pair_projector + ...
    counts.source_scatter + counts.u_cross + counts.u_peanut;

sparse_cols = 0;
sparse_cols = sparse_cols + double(plan.direct_source_corr)*n_source_cols;
sparse_cols = sparse_cols + double(plan.need_u_corr_matrix)*n_u_cols;
sparse_cols = sparse_cols + double(plan.need_pair_nonp_matrix)*n_source_cols;
sparse_cols = sparse_cols + double(plan.combined_pair_u_corr)*n_source_cols;
sparse_cols = sparse_cols + double(plan.need_pair_projector)*pair_total_rows;
sparse_cols = sparse_cols + double(plan.need_source_scatter)*(2*pair_total_rows);
sparse_cols = sparse_cols + double(plan.need_u_cross)*n_u_cols;
sparse_cols = sparse_cols + double(plan.need_u_peanut)*pair_total_rows;

auxiliary_bytes = 8*counts.source_scatter_indices + ...
    8*counts.dense_pair_projector;

estimate = struct();
estimate.matrix_plan = plan;
estimate.counts = counts;
estimate.sparse_cols = sparse_cols;
estimate.estimated_sparse_bytes = 16*counts.total + 8*sparse_cols;
estimate.estimated_auxiliary_bytes = auxiliary_bytes;
estimate.estimated_build_bytes = 24*counts.total;
estimate.estimated_peak_bytes = estimate.estimated_sparse_bytes + ...
    estimate.estimated_build_bytes + estimate.estimated_auxiliary_bytes;
estimate.estimated_sparse_MB = estimate.estimated_sparse_bytes/1024^2;
estimate.estimated_auxiliary_MB = estimate.estimated_auxiliary_bytes/1024^2;
estimate.estimated_build_MB = estimate.estimated_build_bytes/1024^2;
estimate.estimated_peak_MB = estimate.estimated_peak_bytes/1024^2;
end
