function stats = initLaplaceBigSparseSolveStats(requested,n_pairs)
%INITLAPLACEBIGSPARSESOLVESTATS Default stats for Laplace big-sparse solve.

stats = struct();
stats.requested = logical(requested);
stats.active = false;
stats.backend = 'global_block_sparse';
stats.reason = 'not_requested';
stats.build_mode = '';
stats.n_pairs = n_pairs;
stats.N_c = 0;
stats.N_check = 0;
stats.used_pair_cache = false;
stats.rotations_used = false;
stats.source_correction = 'factored_structured';
stats.field_correction = 'factored_sparse';
stats.local_pair_nonp_entries = 0;
stats.local_u_entries = 0;
stats.local_qv_entries = 0;
stats.nnz_u = 0;
stats.nnz_u_cross = 0;
stats.nnz_u_peanut = 0;
stats.nnz_u_qv = 0;
stats.nnz_qv = 0;
stats.nnz_pair_nonp = 0;
stats.nnz_source_scatter = 0;
stats.nnz_qv_scatter = 0;
stats.estimated_sparse_bytes = 0;
stats.estimated_auxiliary_bytes = 0;
stats.estimated_build_bytes = 0;
stats.estimated_peak_bytes = 0;
stats.estimated_sparse_MB = 0;
stats.estimated_auxiliary_MB = 0;
stats.estimated_build_MB = 0;
stats.estimated_peak_MB = 0;
stats.build_time = 0;
if requested
    stats.reason = 'not_prepared';
end
end
