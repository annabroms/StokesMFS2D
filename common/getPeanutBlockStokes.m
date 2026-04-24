function [DC, Y] = getPeanutBlockStokes(rin_pair_c, rin_pair_f, ...
    rout_peanut, Lc_pair, pair_moment_map, pair_rbm_map, ...
    pair_moment_gram, svd_opts)
%GETPEANUTBLOCKSTOKES Compute the Stokes peanut compression factors.
%
% Syntax:
%   [DC, Y] = getPeanutBlockStokes(rin_pair_c,rin_pair_f,rout_peanut, ...
%       Lc_pair,pair_moment_map,pair_rbm_map,pair_moment_gram)
%   [DC, Y] = getPeanutBlockStokes(rin_pair_c,rin_pair_f,rout_peanut, ...
%       Lc_pair,pair_moment_map,pair_rbm_map,pair_moment_gram,svd_opts)
%
% Inputs:
%   rin_pair_c        - coarse pair proxy sources
%   rin_pair_f        - fine pair proxy sources
%   rout_peanut       - peanut collocation targets
%   Lc_pair           - coarse pair projector for the mobility solve
%   pair_moment_map   - 6-by-(2*nsrc) pair force/torque moment map
%   pair_rbm_map      - (2*nsrc)-by-6 pair rigid-body source map
%   pair_moment_gram  - 6-by-6 Gram matrix pair_moment_map*pair_rbm_map
%
% Outputs:
%   DC - Product of the coarse peanut left singular vectors with the
%        projected fine-grid evaluation matrix.
%   Y  - Right pseudoinverse factor returned by getPseudoFactors.
%
% Anna Broms, Feb 12, 2026

if nargin < 8 || isempty(svd_opts)
    svd_opts = struct();
end

mu = 1;

Nf = stokSLPmat(rin_pair_f,rout_peanut,mu);
project = ~isempty(pair_moment_map);
if project
    validate_projection_inputs(Lc_pair,pair_moment_map,pair_rbm_map, ...
        pair_moment_gram);
    Ntot = Nf - (Nf*pair_rbm_map)/pair_moment_gram * pair_moment_map;
else
    Ntot = Nf;
end

Npeanut = stokSLPmat(rin_pair_c,rout_peanut,mu);
if project
    Npeanut = Npeanut*Lc_pair;
end

tol = 1e-14;
svd_opts_local = svd_opts;
if logical(getOptField(svd_opts,'left_weight',false))
    svd_opts_local.row_weights = getPeriodicCurveWeights(rout_peanut);
end

[Y,U] = getPseudoFactors(Npeanut,tol,0,svd_opts_local);
DC = U'*Ntot;
end

function validate_projection_inputs(Lc_pair,pair_moment_map,pair_rbm_map, ...
    pair_moment_gram)
if isempty(Lc_pair) || isempty(pair_moment_map) || isempty(pair_rbm_map) || ...
        isempty(pair_moment_gram)
    error('getPeanutBlockStokes:MissingProjectionData', ...
        ['Projected peanut compression requires Lc_pair, ', ...
         'pair_moment_map, pair_rbm_map, and pair_moment_gram.']);
end
end
