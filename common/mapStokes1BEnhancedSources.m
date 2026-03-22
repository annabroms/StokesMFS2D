function [lambda_x_raw,lambda_y_raw,lambda_body,lambda_px,lambda_py] = mapStokes1BEnhancedSources(tau,geom,basis)
%MAPSTOKES1BENHANCEDSOURCES Map bodywise boundary data to enhanced Stokeslet sources.
%
% Syntax:
%   [lambda_x_raw,lambda_y_raw,lambda_body,lambda_px,lambda_py] = ...
%       mapStokes1BEnhancedSources(tau,geom,basis)
%
% Inputs:
%   tau   - Stacked boundary data [tau_x; tau_y] in bodywise ordering.
%   geom  - Geometry struct from prepareStokes1BEnhanced.
%   basis - Basis struct from prepareStokes1BEnhanced.
%
% Outputs:
%   
%   lambda_x_raw,
%   lambda_y_raw       - Stacked non-projected source strengths
%                        in bodywise ordering.
%   lambda_body        - Non-projected local stacked [x; y] source vector
%                        for each body.
%   lambda_px, lambda_py - Projected stacked Stokeslet source strengths in
%                        bodywise ordering. These are needed for mobility only.

tau = tau(:);
P = numel(geom.q);

[~,lambda_x_raw,lambda_y_raw,lambda_body] = ...
    applyStokes1BEnhancedPseudo(tau,geom,basis);

if geom.solve_resistance
    lambda_px = [];
    lambda_py = [];
else
    lambda_px = zeros(geom.total_source_count,1);
    lambda_py = zeros(geom.total_source_count,1);
end

for k = 1:P
    if ~geom.solve_resistance
        lambda_k = lambda_body{k};
        ns = geom.source_count(k);
        sx = geom.source_x_ind{k};
        sy = geom.source_y_ind{k};
        lambda_k_proj = projectOutRigid2D(lambda_k,geom.rin_body{k},geom.q(k));
        lambda_px(sx) = lambda_k_proj(1:ns);
        lambda_py(sy) = lambda_k_proj(ns+1:end);
    end
end

end
