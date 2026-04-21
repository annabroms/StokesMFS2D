function plan = getMobPeanutBigSparseMatrixPlan(opt)
%GETMOBPEANUTBIGSPARSEMATRIXPLAN Select matrices for big-sparse matvec.
%
% Source correction:
%   direct  -> M_source_corr
%   factored sparse -> M_pair_nonp, M_pair_projector, M_source_scatter
%   factored structured -> M_pair_nonp, P_pair, source_scatter_rows
%
% Velocity correction:
%   direct  -> M_u_corr
%   factored sparse -> M_pair_nonp, M_pair_projector, M_u_cross, M_u_peanut
%   factored structured -> M_pair_nonp, P_pair, M_u_cross, M_u_peanut
%
% Optional comparison mode:
%   combined pair/u -> M_pair_nonp_u_corr = [M_pair_nonp; M_u_corr]

if nargin < 1 || isempty(opt)
    opt = struct();
end

plan = struct();
plan.direct_source_corr = logical(getOptField(opt, ...
    'big_sparse_direct_source_corr',false));
plan.direct_u_corr = logical(getOptField(opt, ...
    'big_sparse_direct_u_corr',false));
plan.structured_apply = logical(getOptField(opt, ...
    'big_sparse_structured_apply',false));

plan.need_factored_source_corr = ~plan.direct_source_corr;
plan.need_factored_u_corr = ~plan.direct_u_corr;
plan.need_pair_nonp = plan.need_factored_source_corr || ...
    plan.need_factored_u_corr;
plan.need_dense_pair_projector = plan.structured_apply && ...
    plan.need_pair_nonp;
plan.need_pair_projector = plan.need_pair_nonp && ...
    ~plan.need_dense_pair_projector;
plan.need_source_scatter_indices = plan.structured_apply && ...
    plan.need_factored_source_corr;
plan.need_source_scatter = plan.need_factored_source_corr && ...
    ~plan.need_source_scatter_indices;
plan.need_u_cross = plan.need_factored_u_corr;
plan.need_u_peanut = plan.need_factored_u_corr;
plan.combined_pair_u_corr = logical(getOptField(opt, ...
    'big_sparse_combine_pair_u_corr',false)) && ...
    plan.need_pair_nonp && plan.direct_u_corr;
plan.need_pair_nonp_matrix = plan.need_pair_nonp && ...
    ~plan.combined_pair_u_corr;
plan.need_u_corr_matrix = plan.direct_u_corr && ...
    ~plan.combined_pair_u_corr;
end
