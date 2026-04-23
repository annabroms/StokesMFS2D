function mode = resolveResBigSparseUCorrMode(opt)
%RESOLVERESBIGSPARSEUCORRMODE Select resistance sparse velocity map layout.

if nargin < 1 || isempty(opt)
    opt = struct();
end

if isfield(opt,'res_big_sparse_u_corr_mode')
    mode = lower(char(opt.res_big_sparse_u_corr_mode));
else
    mode = 'factored';
end

switch mode
    case {'factored','combined'}
        % supported
    case {'one_matrix','combined_sparse'}
        mode = 'combined';
    otherwise
        error('resolveResBigSparseUCorrMode:BadMode', ...
            'opt.res_big_sparse_u_corr_mode must be ''factored'' or ''combined''.');
end
end
