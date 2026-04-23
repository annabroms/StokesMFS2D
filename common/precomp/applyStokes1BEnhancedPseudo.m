function [lambda,lambda_x,lambda_y,lambda_body] = applyStokes1BEnhancedPseudo(u,geom,basis)
%APPLYSTOKES1BENHANCEDPSEUDO Apply the bodywise enhanced 1B pseudoinverse.
%
% Syntax:
%   [lambda,lambda_x,lambda_y,lambda_body] = ...
%       applyStokes1BEnhancedPseudo(u,geom,basis)
%
% Inputs:
%   u     - Stacked target data [u_x; u_y] in bodywise ordering.
%   geom  - Geometry struct from prepareStokes1BEnhanced.
%   basis - Basis struct from prepareStokes1BEnhanced.
%
% Outputs:
%   lambda      - Stacked source strengths [lambda_x; lambda_y].
%   lambda_x    - x-strengths in bodywise source ordering.
%   lambda_y    - y-strengths in bodywise source ordering.
%   lambda_body - Cell array with bodywise stacked [x; y] strengths.
%
% Notes:
%   Each body block is applied in the factorized order
%
%       Aii^+ * u_i = Y_i * (U_i * u_i),
%
%   where prepareStokes1BEnhanced stores U_i as the transposed left-singular
%   factor so the multiplication order stays explicit and shared.

u = u(:);
expected_len = 2*geom.total_target_count;
if numel(u) ~= expected_len
    error('applyStokes1BEnhancedPseudo:badInputSize', ...
        'u must have length 2*geom.total_target_count.');
end

P = numel(geom.q);
lambda_x = zeros(geom.total_source_count,1);
lambda_y = zeros(geom.total_source_count,1);
lambda_body = cell(P,1);

for k = 1:P
    uk = [u(geom.target_x_ind{k}); u(geom.target_y_ind{k})];
    lambda_k = basis.Y{k}*(basis.U{k}*uk);
    lambda_body{k} = lambda_k;

    ns = geom.source_count(k);
    sx = geom.source_x_ind{k};
    lambda_x(sx) = lambda_k(1:ns);
    lambda_y(sx) = lambda_k(ns+1:end);
end

lambda = [lambda_x; lambda_y];

end
