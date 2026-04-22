function estimate = estimateMobPeanutBigSparseRamStokes(P,N_c,N_check,n_pairs,opt)
%ESTIMATEMOBPEANUTBIGSPARSERAMSTOKES Estimate big-sparse matrix memory.
%
% estimate = estimateMobPeanutBigSparseRamStokes(P,N_c,N_check,n_pairs,opt)
% returns conservative byte counts for the sparse matrices and temporary
% sparse-entry arrays used by buildMobPeanutBigSparseStokes. Auxiliary
% arrays cover the dense local pair projector and indexed source scatter
% rows used by the fixed source correction. The explicit
% retained-pair-basis estimate distinguishes precomputed mode, which keeps
% dense Cmap/Cmap_FU pair maps, from streaming mode, which discards them.
%
% opt.big_sparse_direct_u_corr selects the velocity correction matrices:
%   true  (default): M_u_corr maps projected one-body sources to u_corr.
%   false: M_u_cross and M_u_peanut keep the velocity correction factored.
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

direct_u_corr = logical(getOptField(opt,'big_sparse_direct_u_corr',true));
build_mode = lower(char(getOptField(opt,'mob_big_sparse_build_mode', ...
    'precomputed')));
switch build_mode
    case {'precomputed','streaming'}
        % supported
    otherwise
        error('estimateMobPeanutBigSparseRamStokes:BadMode', ...
            ['opt.mob_big_sparse_build_mode must be ''precomputed'' ', ...
             'or ''streaming''.']);
end
chunk_pairs = max(1,round(getOptField(opt,'mob_big_sparse_chunk_pairs',8)));
chunk_pairs = min(max(1,n_pairs),chunk_pairs);
if direct_u_corr
    velocity_correction = 'direct_sparse';
else
    velocity_correction = 'factored_sparse';
end
if logical(getOptField(opt,'use_matrix_free_Lc_pair',true))
    projector_mode = 'matrix_free_equivalent';
else
    projector_mode = 'dense_Lc_pair';
end
plan = struct('source_correction','factored_structured', ...
    'velocity_correction',velocity_correction, ...
    'direct_u_corr',direct_u_corr, ...
    'build_mode',build_mode, ...
    'chunk_pairs',chunk_pairs, ...
    'projector_mode',projector_mode);

n_coarse = P*N_c;
n_source_cols = 2*n_coarse;
n_u_cols = 2*n_coarse;
pair_rows = 4*N_c;
pair_total_rows = n_pairs*pair_rows;

u_entries = n_pairs*(4*N_check)*(4*N_c);
pair_nonp_entries = n_pairs*pair_rows*(4*N_c);
rbm_entries = n_pairs*6*(4*N_c);
source_scatter_index_entries = n_pairs*(8*N_c);
dense_pair_projector_entries = pair_rows*pair_rows;

counts = struct();
counts.u = double(direct_u_corr)*u_entries;
counts.pair_nonp = pair_nonp_entries;
counts.rbm = rbm_entries;
counts.source_scatter_indices = source_scatter_index_entries;
counts.dense_pair_projector = dense_pair_projector_entries;
counts.u_cross = double(~direct_u_corr)*u_entries;
counts.u_peanut = double(~direct_u_corr)*u_entries;
counts.total = counts.u + counts.pair_nonp + counts.u_cross + ...
    counts.u_peanut + counts.rbm;
if strcmp(build_mode,'precomputed')
    counts.retained_pair_basis = pair_nonp_entries + rbm_entries;
else
    counts.retained_pair_basis = 0;
end

sparse_cols = n_source_cols + double(direct_u_corr)*n_u_cols + ...
    double(~direct_u_corr)*(n_u_cols + pair_total_rows) + ...
    n_source_cols;
local_u_entries = (4*N_check)*(4*N_c);
local_pair_nonp_entries = pair_rows*(4*N_c);
local_rbm_entries = 6*(4*N_c);
chunk_entries = chunk_pairs*(local_pair_nonp_entries + ...
    double(direct_u_corr)*local_u_entries + ...
    double(~direct_u_corr)*2*local_u_entries + local_rbm_entries);

auxiliary_bytes = 8*counts.source_scatter_indices + ...
    8*counts.dense_pair_projector;
retained_pair_basis_bytes = 8*counts.retained_pair_basis;

estimate = struct();
estimate.matrix_plan = plan;
estimate.counts = counts;
estimate.sparse_cols = sparse_cols;

estimate.big_sparse_matrix_bytes = 16*counts.total + 8*sparse_cols;
estimate.big_sparse_auxiliary_bytes = auxiliary_bytes;
estimate.big_sparse_build_bytes = 24*chunk_entries;
estimate.big_sparse_peak_bytes = estimate.big_sparse_matrix_bytes + ...
    estimate.big_sparse_build_bytes + estimate.big_sparse_auxiliary_bytes;
estimate.retained_pair_basis_bytes = retained_pair_basis_bytes;
estimate.solver_precompute_peak_bytes = estimate.big_sparse_peak_bytes + ...
    estimate.retained_pair_basis_bytes;
end
