function estimate = estimateResPeanutBigSparseRamStokes(P,N_c,N_check,n_pairs,opt)
%ESTIMATERESPEANUTBIGSPARSERAMSTOKES Estimate resistance sparse-map memory.

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
combined_u_corr = strcmp(u_corr_mode,'combined');
velocity_correction = [u_corr_mode '_sparse'];

n_coarse = P*N_c;
pair_rows = 4*N_c;
pair_total_rows = n_pairs*pair_rows;

local_pair = pair_rows*(4*N_c);
local_u = (4*N_check)*(4*N_c);
local_ft = 6*(4*N_c);

counts = struct();
counts.pair = n_pairs*local_pair;
counts.u = double(combined_u_corr)*n_pairs*2*local_u;
counts.u_cross = double(strcmp(u_corr_mode,'factored'))*n_pairs*local_u;
counts.u_peanut = double(strcmp(u_corr_mode,'factored'))*n_pairs*local_u;
counts.ft = n_pairs*local_ft;
counts.source_scatter_indices = n_pairs*(4*N_c);
counts.velocity_scatter_indices = n_pairs*(4*N_check);
counts.ft_scatter_indices = n_pairs*6;
counts.total = counts.pair + counts.u + counts.u_cross + counts.u_peanut + ...
    counts.ft;

n_source_cols = 2*n_coarse;
sparse_cols = n_source_cols + n_source_cols + n_source_cols + ...
    pair_total_rows;
chunk_entries = chunk_pairs*(local_pair + 2*local_u + local_ft);

estimate = struct();
estimate.matrix_plan = struct( ...
    'source_correction','direct_sparse_scatter', ...
    'velocity_correction',velocity_correction, ...
    'u_corr_mode',u_corr_mode, ...
    'combined_u_corr',combined_u_corr, ...
    'force_torque_correction','direct_sparse', ...
    'chunk_pairs',chunk_pairs);
estimate.counts = counts;
estimate.big_sparse_matrix_bytes = 16*counts.total + 8*sparse_cols;
estimate.big_sparse_auxiliary_bytes = 8*(counts.source_scatter_indices + ...
    counts.velocity_scatter_indices + counts.ft_scatter_indices);
estimate.big_sparse_build_bytes = 24*chunk_entries;
estimate.big_sparse_peak_bytes = estimate.big_sparse_matrix_bytes + ...
    estimate.big_sparse_auxiliary_bytes + estimate.big_sparse_build_bytes;
end
