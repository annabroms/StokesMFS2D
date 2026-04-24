function [lam_c_x, lam_c_y, lambda_self_blocks, u_corr] = ...
    transform_res_peanut_big_sparse_stokes(tau,geom,basis)
%TRANSFORM_RES_PEANUT_BIG_SPARSE_STOKES Apply one-body pseduoinverse plus 
% solve-grid sparse pair maps.

if ~isfield(basis,'big_sparse') || isempty(basis.big_sparse)
    error('transform_res_peanut_big_sparse_stokes:MissingBigSparse', ...
        'basis.big_sparse is required for opt.use_big_sparse=1.');
end
if ~isequal(geom.rcheck,geom.rvec_out)
    error('transform_res_peanut_big_sparse_stokes:UnsupportedGrid', ...
        'The resistance big sparse transform is only supported on the solve grid.');
end

opt = geom.opt;
rvec_out = geom.rvec_out;

P = opt.P;
N_c = opt.N_c;
N_large = numel(rvec_out)/P;
PM = numel(rvec_out);
N_check = PM/P;
n_coarse = P*N_c;

%% For the one-body basis:
% Map boundary data through the one-body pseudoinverse exactly as in the
% loop transform but batch them as all one-body blocks are the same. 
% This is faster than a loop over P, see
% experiments/apr23_resistance_one_body_batch_compare.m
% These sources also drive all pair corrections.

tau_x = reshape(tau(1:PM),N_large,P);
tau_y = reshape(tau(PM+1:2*PM),N_large,P);
tau_blocks = [tau_x; tau_y];

step_blocks = basis.U{1}*tau_blocks;
lambda_self_blocks = basis.Y{1}*step_blocks;

lam_self_x = reshape(lambda_self_blocks(1:N_c,:),[],1);
lam_self_y = reshape(lambda_self_blocks(N_c+1:2*N_c,:),[],1);

lambda_self = [lam_self_x; lam_self_y];

%% Apply all pair corrections through a fixed sparse map.
pair_tau = basis.big_sparse.M_pair*lambda_self;
if isfield(basis.big_sparse,'M_source_corr') && ...
        ~isempty(basis.big_sparse.M_source_corr)
    lambda_c = lambda_self + basis.big_sparse.M_source_corr*lambda_self;
else
    lambda_c = addSourceCorrectionsInPairOrder(lambda_self, ...
        basis.big_sparse.source_scatter_rows,pair_tau,4*N_c);
end
lam_c_x = lambda_c(1:n_coarse);
lam_c_y = lambda_c(n_coarse+1:2*n_coarse);

% Correct the velocity field for all pairs. Combined mode uses one
% sparse matrix on [lambda_self; pair_tau].
switch getVelocityCorrectionMode(opt,basis.big_sparse)
    case 'combined'
        if ~isfield(basis.big_sparse,'M_u_corr')
            error('transform_res_peanut_big_sparse_stokes:MissingUCorr', ...
                'Combined resistance velocity correction requires M_u_corr.');
        end
        pair_u_corr = basis.big_sparse.M_u_corr*[lambda_self; pair_tau];
    otherwise
        if ~isfield(basis.big_sparse,'M_u_cross') || ...
                ~isfield(basis.big_sparse,'M_u_peanut')
            error('transform_res_peanut_big_sparse_stokes:MissingFactoredUCorr', ...
                ['Factored resistance velocity correction requires ', ...
                 'M_u_cross and M_u_peanut.']);
        end
        pair_u_corr = basis.big_sparse.M_u_cross*lambda_self - ...
            basis.big_sparse.M_u_peanut*pair_tau;
end
    u_corr = addVectorCorrectionsInPairOrder(zeros(2*numel(geom.rcheck),1), ...
        basis.big_sparse.velocity_scatter_rows,pair_u_corr,4*N_check);

end

function mode = getVelocityCorrectionMode(opt,big_sparse)
    if isfield(big_sparse,'matrix_plan') && ...
            isfield(big_sparse.matrix_plan,'u_corr_mode')
        mode = big_sparse.matrix_plan.u_corr_mode;
    else
        mode = resolveResBigSparseUCorrMode(opt);
    end
end

function lambda_c = addSourceCorrectionsInPairOrder(lambda_self,rows, ...
    pair_tau,pair_rows)
    n_pairs = numel(pair_tau)/pair_rows;
    lambda_c = lambda_self;
    for row = 1:n_pairs
        idx = (row-1)*pair_rows+1:row*pair_rows;
        lambda_c(rows(idx)) = lambda_c(rows(idx)) + pair_tau(idx);
    end
end

function out = addVectorCorrectionsInPairOrder(out,rows,pair_values,pair_rows)
    n_pairs = numel(pair_values)/pair_rows;
    for row = 1:n_pairs
        idx = (row-1)*pair_rows+1:row*pair_rows;
        out(rows(idx)) = out(rows(idx)) + pair_values(idx);
    end
end
