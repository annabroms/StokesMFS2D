function [Y, U] = getPseudoFactors(N, tol, visualise, svd_opts)
%GETPSEUDOFACTORS Computes factors that provide the matrix psuedoinverse from 
% a truncated SVD 
%
% Syntax:
%   [Y, U] = getPseudoFactors(N, tol, visualise)
%   [Y, U] = getPseudoFactors(N, tol, visualise, svd_opts)
%
% Inputs:
%   N         - Input matrix for which to compute a pseudoinverse
%   tol       - Relative truncation threshold; singular values σ are kept if σ > max(σ) * tol
%   visualise - Logical flag: plot singular values (true/false)
%   svd_opts  - Optional struct. Supported fields:
%               column_weight : scale each column by the reciprocal of its
%                               2-norm before the SVD, while returning
%                               factors for the original variable.
%               left_weight   : scale each row by sqrt(row_weights)
%                               before the SVD, while returning factors
%                               that still apply as x = Y*(U'*rhs).
%               row_weights   : positive row weights used when
%                               left_weight = true. Can have length
%                               size(N,1) or size(N,1)/2, in which case
%                               the weights are repeated for x/y Stokes
%                               row blocks.
%
% Outputs:
%   U - Matrix of left singular vectors corresponding to retained singular values
%   Y - Product VS⁺, where:
%         - S⁺ is a diagonal matrix with entries 1/σ for retained singular values
%         - V contains the corresponding right singular vectors
%
% Description:
%   Computes a truncated SVD of matrix N and returns factor matrices U and Y such that:
%       N⁺ ≈ Y * U'
%   This allows efficient and backward-stable application of the pseudoinverse without explicitly forming it.
%
% Notes:
%   - Only singular values greater than max(σ) * tol are retained
%   - Intended for use in stable pseudoinverse application (e.g., solving least-squares problems)
%
% Anna Broms 4 April 2025

if nargin == 0
    test_getPseudoFactors_weighting;
    return
end

% - With N the discretized Stokeslets, 
%with tol = eps, we throw away only the last singval, corresponding to the
%non-trivial null-space (normal direction for single layer...). Does that make sense?
%yes! testing ra = size(N,2) leads to a horribly irregular 1-body basis. 

%econ = 0; %What does econ do here?? 

% if econ
%     [UU,S,V] = svd(N,"econ");
%     S = diag(S);
%     iS = 1./S; % rank
%     Y = V*diag(iS); 
%     Bi1 = UU; 
%else

if nargin < 3
    visualise = 0; 
end
if nargin < 4 || isempty(svd_opts)
    svd_opts = struct();
end

column_weight = logical(getOptField(svd_opts,'column_weight',false));
left_weight = logical(getOptField(svd_opts,'left_weight',false));

if column_weight || left_weight
    N_svd = N;
    if column_weight
        col_norms = vecnorm(N,2,1);
        col_scale = ones(size(col_norms));
        nz = col_norms > 0;
        col_scale(nz) = 1./col_norms(nz);
        N_svd = bsxfun(@times,N_svd,col_scale);
    else
        col_scale = [];
    end

    if left_weight
        row_weights = getOptField(svd_opts,'row_weights',[]);
        row_scale = get_row_scale(row_weights,size(N,1));
        N_svd = bsxfun(@times,row_scale,N_svd);
    else
        row_scale = [];
    end
else
    col_scale = [];
    row_scale = [];
    N_svd = N;
end

[UU,S,V] = svd(N_svd);
S = diag(S);

%use relative tolerance 
if isempty(S)
    ra = 0;
else
    ra = sum(S>max(S)*tol);
end


if visualise
   
    figure(57)
   % clf;
    semilogy(S,'o-');
    hold on
    semilogy(ra*ones(1,2),logspace(-15,5,2),'r--')
    if isempty(S) || min(S) == 0
        c = inf;
    else
        c = max(S)/min(S);  %Condition number
    end
    if isempty(S)
        smax = 0;
    else
        smax = max(S);
    end
    str = sprintf('Self condition number %1.3e, max sing %1.3e',c,smax);
    title(str,'interpreter','latex');
    grid on
    xlabel('$j$','interpreter','latex')
    ylabel('$\sigma_j$','interpreter','latex')

    % figure(56)
    % SV = SS*V';
    % semilogy(abs(SV(:,end-50:end)'))
    % 
    % 
    % figure(56)
    % semilogy(abs(UU(end-3:end,:)'))

    % figure(57)
    % semilogy(abs(diff(S)));
    % hold on
    % title('Decay rate of sing vals','interpreter','latex')
    % 
    % 
    % figure(58)
    % semilogy(abs(diff(S)./S(2:end)));
    % hold on
    % title('Relative decay rate of sing vals','interpreter','latex')
end

S = S(1:ra);  %get pseudoinverse of S
iS = 1./S; 
Y = V(:,1:ra)*diag(iS);
if column_weight
    Y = bsxfun(@times,col_scale(:),Y);
end
U = UU(:,1:ra);
if left_weight
    U = bsxfun(@times,row_scale,U);
end





end

function row_scale = get_row_scale(row_weights,nrows)
row_weights = row_weights(:);
if isempty(row_weights)
    error('getPseudoFactors:MissingRowWeights', ...
        'svd_opts.row_weights must be supplied when left_weight = true.');
end

if numel(row_weights) == nrows
    w_full = row_weights;
elseif mod(nrows,2) == 0 && numel(row_weights) == nrows/2
    w_full = repmat(row_weights,2,1);
else
    error('getPseudoFactors:InvalidRowWeights', ...
        ['svd_opts.row_weights must have length size(N,1) or size(N,1)/2 ', ...
         'for Stokes x/y row blocks.']);
end

if any(~isfinite(w_full))
    error('getPseudoFactors:InvalidRowWeights', ...
        'svd_opts.row_weights must be finite.');
end

positive = w_full > 0;
if ~all(positive)
    if any(positive)
        floor_w = eps(max(1,max(w_full(positive))));
    else
        floor_w = 1;
    end
    w_full(~positive) = floor_w;
end

row_scale = sqrt(w_full);
end

function test_getPseudoFactors_weighting
tol = 1e-12;
visualise = 0;

cases = cell(2,1);
cases{1} = [1e-10 1 5e3; -2e-10 2 -1e3; 3e-10 -4 2e3; 4e-10 0.5 7e2];
cases{2} = [1 0 3; 2 0 -1; -1 0 4; 0.5 0 2];
row_weights = {
    [1e-3; 3; 0.2; 7], ...
    [4; 1e-2; 5; 0.5]
};
b_list = {[1; -2; 3; -1], [0.3; 0.8; -1.1; 2.0]};
modes = { ...
    struct('name','column','column_weight',true,'left_weight',false), ...
    struct('name','left','column_weight',false,'left_weight',true), ...
    struct('name','both','column_weight',true,'left_weight',true)};

for icase = 1:numel(cases)
    N = cases{icase};
    col_norms = vecnorm(N,2,1);
    col_scale = ones(size(col_norms));
    nz = col_norms > 0;
    col_scale(nz) = 1./col_norms(nz);
    D = diag(col_scale);
    W = diag(sqrt(row_weights{icase}));

    for imode = 1:numel(modes)
        svd_opts = modes{imode};
        svd_opts.row_weights = row_weights{icase};

        [Y,U] = getPseudoFactors(N,tol,visualise,svd_opts);
        assert(all(isfinite(Y(:))),'Weighted Y contains non-finite entries.');
        assert(all(isfinite(U(:))),'Weighted U contains non-finite entries.');

        if svd_opts.column_weight
            Dmode = D;
        else
            Dmode = eye(size(N,2));
        end
        if svd_opts.left_weight
            Wmode = W;
        else
            Wmode = eye(size(N,1));
        end

        Nwd = Wmode*N*Dmode;
        sing = svd(Nwd);
        if isempty(sing)
            pinv_tol = 0;
        else
            pinv_tol = max(sing)*tol;
        end

        for ib = 1:numel(b_list)
            b = b_list{ib};
            x_fac = Y*(U'*b);
            x_ref = Dmode*pinv(Nwd,pinv_tol)*Wmode*b;
            err = norm(x_fac-x_ref);
            assert(err <= 1e-11*max(1,norm(x_ref)), ...
                'Weighted pseudoinverse action mismatch in case %d, mode %s.', ...
                icase,svd_opts.name);
            assert(all(isfinite(x_fac)),'Weighted solve produced non-finite entries.');
        end
    end
end

fprintf('getPseudoFactors weighting self-test passed.\n');
end
