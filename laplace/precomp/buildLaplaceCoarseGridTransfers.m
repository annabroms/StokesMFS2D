function transfers = buildLaplaceCoarseGridTransfers( ...
    rbase_in_c,rbase_in_cmap,rotations,R,a_c,tol,svd_opts)
%BUILDLAPLACECOARSEGRIDTRANSFERS Rectangular, field-preserving source maps.
%
% Syntax:
%   transfers = buildLaplaceCoarseGridTransfers( ...
%       rbase_in_c,rbase_in_cmap,rotations,R,a_c,tol,svd_opts)
%
% The solver coefficients live on rbase_in_c (N_c nodes), whereas a
% canonical pair Cmap may use rbase_in_cmap (N_cmap nodes).  For every
% pair-axis rotation rotations(k), this routine constructs
%
%   to_cmap{k}   : N_cmap-by-N_c, actual grid -> canonical Cmap grid,
%   from_cmap{k} : N_c-by-N_cmap, canonical Cmap grid -> actual grid.
%
% Both maps refit the represented one-body Laplace field on an oversampled
% particle boundary.  They preserve the sum of source strengths exactly;
% this is the logarithmic far-field coefficient and hence the body charge.
% The destination pseudoinverses are factored only once, while the forward
% evaluation matrices depend on the pair angle.

if nargin < 4 || isempty(R)
    R = 1;
end
if nargin < 5 || isempty(a_c)
    a_c = 1.2;
end
if nargin < 6 || isempty(tol)
    tol = 1e-14;
end
if nargin < 7 || isempty(svd_opts)
    svd_opts = struct();
end

rbase_in_c = rbase_in_c(:);
rbase_in_cmap = rbase_in_cmap(:);
rotations = rotations(:);
if isempty(rotations)
    rotations = 1;
end

N_c = numel(rbase_in_c);
N_cmap = numel(rbase_in_cmap);
validateattributes(N_c,{'numeric'},{'scalar','integer','positive'});
validateattributes(N_cmap,{'numeric'},{'scalar','integer','positive'});
validateattributes(R,{'numeric'},{'scalar','real','positive','finite'});
validateattributes(a_c,{'numeric'},{'scalar','real','>=',1,'finite'});
validateattributes(tol,{'numeric'},{'scalar','real','positive','finite'});

% Use enough targets to resolve the larger destination space.  A half-step
% angular offset keeps validation grids in callers naturally disjoint from
% this training grid.
N_train = max(ceil(a_c*max(N_c,N_cmap)),max(N_c,N_cmap)+1);
t_train = (0:N_train-1)'*(2*pi/N_train) + pi/N_train;
rout_train = R*exp(1i*t_train);

context_c = prepare_destination(rbase_in_c,rout_train,tol,svd_opts);
context_cmap = prepare_destination(rbase_in_cmap,rout_train,tol,svd_opts);

nrot = numel(rotations);
to_cmap = cell(nrot,1);
from_cmap = cell(nrot,1);
for k = 1:nrot
    rot = rotations(k);
    if abs(rot) == 0
        rot = 1;
    else
        rot = rot/abs(rot);
    end

    % Move the actual solver grid into the canonical pair frame.
    A_from = lapSLPmat(conj(rot)*rbase_in_c,rout_train);
    to_cmap{k} = finish_transfer(context_cmap,A_from,N_c);

    % Rotate the canonical Cmap grid into the actual physical frame.
    A_from = lapSLPmat(rot*rbase_in_cmap,rout_train);
    from_cmap{k} = finish_transfer(context_c,A_from,N_cmap);
end

transfers = struct();
transfers.enabled = true;
transfers.method = 'charge_preserving_field_refit';
transfers.N_c = N_c;
transfers.N_cmap = N_cmap;
transfers.N_train = N_train;
transfers.tol = tol;
transfers.to_cmap = to_cmap;
transfers.from_cmap = from_cmap;
end

function context = prepare_destination(rsrc,rout,tol,svd_opts)
n = numel(rsrc);
A = lapSLPmat(rsrc,rout);
P0 = eye(n) - ones(n)/n;
[Y,U] = getPseudoFactors(A*P0,tol,0,svd_opts);
context = struct('A',A,'P0',P0,'Y',Y,'Ut',U');
end

function F = finish_transfer(context,A_from,n_from)
n_to = size(context.A,2);

% C carries the constant Fourier mode and maps total source strength
% exactly: ones(1,n_to)*C = ones(1,n_from).
C = (ones(n_to,1)/n_to)*ones(1,n_from);
rhs = A_from - context.A*C;
F = C + context.P0*(context.Y*(context.Ut*rhs));

% Remove roundoff in the charge constraint without changing the field fit
% at meaningful precision.
charge_error = ones(1,n_from) - sum(F,1);
F = F + (ones(n_to,1)/n_to)*charge_error;
end
