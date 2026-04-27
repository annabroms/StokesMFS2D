function estimate = estimateMobPeanutBigSparseRam(P,N_c,N_check,n_pairs,opt)
%ESTIMATEMOBPEANUTBIGSPARSERAM Estimate big-sparse matrix memory.
%
% estimate = estimateMobPeanutBigSparseRam(P,N_c,N_check,n_pairs,opt)
% returns conservative byte counts for the sparse matrices and temporary
% sparse-entry arrays used by buildMobPeanutBigSparse. Auxiliary
% arrays cover retained source-correction helpers on the legacy
% mob_sparse_map_coarse=0 path. The explicit
% retained-pair-basis estimate distinguishes precomputed mode, which keeps
% dense pair maps and, when boundary postprocessing needs it, the explicit
% UB/YB pair source factors, from streaming mode, which discards them.
%
% opt.mob_sparse_map_coarse selects the mobility coarse source-correction
% path:
%   false: keep the legacy factored source correction
%                    M_pair_nonp -> P_pair -> source_scatter_rows.
%   true  (default): build a direct global M_source_corr instead.
%
% The legacy shared alias opt.sparse_map_coarse is still accepted.
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
sparse_map_coarse = logical(getOptField(opt,'mob_sparse_map_coarse',...
    getOptField(opt,'sparse_map_coarse',true)));
build_mode = lower(char(getOptField(opt,'mob_big_sparse_build_mode', ...
    'precomputed')));
switch build_mode
    case {'precomputed','streaming'}
        % supported
    otherwise
        error('estimateMobPeanutBigSparseRam:BadMode', ...
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
if sparse_map_coarse
    source_correction = 'direct_sparse_global';
elseif logical(getOptField(opt,'use_matrix_free_Lc_pair',true))
    source_correction = 'factored_structured';
    projector_mode = 'matrix_free_equivalent';
else
    source_correction = 'factored_structured';
    projector_mode = 'dense_Lc_pair';
end
if sparse_map_coarse
    if direct_u_corr
        projector_mode = 'offline_source_projection';
    else
        projector_mode = 'offline_source_and_pair_projection';
    end
end
plan = struct('source_correction','factored_structured', ...
    'sparse_map_coarse',sparse_map_coarse, ...
    'velocity_correction',velocity_correction, ...
    'direct_u_corr',direct_u_corr, ...
    'build_mode',build_mode, ...
    'chunk_pairs',chunk_pairs, ...
    'projector_mode',projector_mode);
plan.source_correction = source_correction;

n_coarse = P*N_c;
n_source_cols = 2*n_coarse;
n_u_cols = 2*n_coarse;
pair_rows = 4*N_c;
pair_total_rows = n_pairs*pair_rows;
N_f = getOptField(opt,'N_f',N_c);
a_f = getOptField(opt,'a_f',1.2);
Nclust = getOptField(opt,'Nclust',0);

local_source = pair_rows*(4*N_c);
u_entries = n_pairs*(4*N_check)*(4*N_c);
pair_nonp_entries = n_pairs*local_source;
source_corr_entries = n_pairs*local_source;
pair_proj_entries = double(sparse_map_coarse && ~direct_u_corr)* ...
    n_pairs*local_source;
rbm_entries = n_pairs*6*(4*N_c);
source_scatter_index_entries = double(~sparse_map_coarse)*n_pairs*(8*N_c);
dense_pair_projector_entries = double(~sparse_map_coarse)*(pair_rows*pair_rows);
full_pair_basis_entries = estimate_full_pair_basis_entries( ...
    n_pairs,N_c,N_f,a_f,Nclust);
retains_full_pair_basis = strcmp(build_mode,'precomputed') && ...
    uses_full_pair_payload(opt);

counts = struct();
counts.u = double(direct_u_corr)*u_entries;
counts.source_corr = double(sparse_map_coarse)*source_corr_entries;
counts.pair_nonp = double(~sparse_map_coarse)*pair_nonp_entries;
counts.pair_proj = pair_proj_entries;
counts.rbm = rbm_entries;
counts.source_scatter_indices = source_scatter_index_entries;
counts.dense_pair_projector = dense_pair_projector_entries;
counts.u_cross = double(~direct_u_corr)*u_entries;
counts.u_peanut = double(~direct_u_corr)*u_entries;
counts.total = counts.u + counts.source_corr + counts.pair_nonp + ...
    counts.pair_proj + counts.u_cross + counts.u_peanut + counts.rbm;
if strcmp(build_mode,'precomputed')
    counts.retained_pair_basis = pair_nonp_entries + rbm_entries + ...
        double(retains_full_pair_basis)*full_pair_basis_entries;
else
    counts.retained_pair_basis = 0;
end

sparse_cols = double(sparse_map_coarse)*n_source_cols + ...
    double(~sparse_map_coarse)*n_source_cols + ...
    double(sparse_map_coarse && ~direct_u_corr)*n_source_cols + ...
    double(direct_u_corr)*n_u_cols + ...
    double(~direct_u_corr)*(n_u_cols + pair_total_rows) + ...
    n_source_cols;
local_u_entries = (4*N_check)*(4*N_c);
local_pair_nonp_entries = pair_rows*(4*N_c);
local_source_corr_entries = local_pair_nonp_entries;
local_pair_proj_entries = double(sparse_map_coarse && ~direct_u_corr)* ...
    local_pair_nonp_entries;
local_rbm_entries = 6*(4*N_c);
chunk_entries = chunk_pairs*( ...
    double(sparse_map_coarse)*local_source_corr_entries + ...
    double(~sparse_map_coarse)*local_pair_nonp_entries + ...
    local_pair_proj_entries + ...
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

function entries = estimate_full_pair_basis_entries(n_pairs,N_c,N_f,a_f,Nclust)
% Estimate the explicit UB/YB pair factors retained by precomputed boundary
% postprocessing. Nclust is a per-side upper bound for enhancement nodes.
if n_pairs == 0
    entries = 0;
    return
end

nout_f = ceil(a_f*N_f);
n_enhance_pair = 2*max(0,Nclust);
n_sources_pair = 2*N_f + n_enhance_pair;
n_targets_pair = 2*nout_f + n_enhance_pair;
rank_bound = min(2*n_sources_pair,2*n_targets_pair);

upf_entries = rank_bound*(4*N_c);
ypf_entries = (2*n_sources_pair)*rank_bound;
entries = n_pairs*(upf_entries + ypf_entries);
end

function tf = uses_full_pair_payload(opt)
N_peanut = getOptField(opt,'N_peanut',0);
use_pair_map = logical(getOptField(opt,'cmap',false));
get_bndry_field = logical(getOptField(opt,'get_bndry_field',false));
self_correct = logical(getOptField(opt,'self_correct',false));
debug = logical(getOptField(opt,'pair_basis_debug',false));
tf = ~((N_peanut > 0) && use_pair_map && self_correct && ...
    ~get_bndry_field && ~debug);
end
