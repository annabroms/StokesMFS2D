function [U,Y,Atot] = getPairBlockStokes(rin_pair,rout_f, ...
    pair_moment_map,pair_rbm_map,pair_moment_gram,pair_target_rbm_map,svd_opts)
%GETPAIRBLOCKSTOKES Build pseudoinverse factors for a 2-body Stokes pair block.
%
% Syntax:
%   [U,Y] = getPairBlockStokes(rin_pair,rout_f,pair_moment_map, ...
%       pair_rbm_map,pair_moment_gram,pair_target_rbm_map)
%   [U,Y] = getPairBlockStokes(rin_pair,rout_f,pair_moment_map, ...
%       pair_rbm_map,pair_moment_gram,pair_target_rbm_map,svd_opts)
%
% Inputs:
%   rin_pair  - complex source locations for the pair fine basis (body 1
%               followed by body 2, including any enhancement points).
%   rout_f    - complex collocation locations for the pair solve
%               (body 1 followed by body 2).
%   pair_moment_map     - 6-by-(2*nsrc) map from pair-ordered source
%                         strengths to per-body force/torque moments.
%   pair_rbm_map        - (2*nsrc)-by-6 map from rigid-body moments back to
%                         pair-ordered source strengths.
%   pair_moment_gram    - 6-by-6 Gram matrix pair_moment_map*pair_rbm_map.
%   pair_target_rbm_map - (2*ntarget)-by-6 rigid-body velocity map on the
%                         pair collocation grid. If empty, no projection is
%                         applied and the raw single-layer operator is used.
%
% Outputs:
%   U, Y  - pseudoinverse factors returned by getPseudoFactors for the
%           assembled pair operator Atot, such that the stabilized solve
%           can be applied as x = Y*(U'*rhs).
%
% Operator assembled in this routine:
%   S    = stokSLPmat(rin_pair, rout_f, 1)
%   Atot = S + (pair_target_rbm_map - (S*pair_rbm_map)/pair_moment_gram) ...
%          * pair_moment_map
%   Atot = S   when the projection data is empty.
%
% See also: getPseudoFactors, stokSLPmat, getKftPair
%
% Anna Broms Feb 12, 2026

if nargin < 7 || isempty(svd_opts)
    svd_opts = struct();
end

mu = 1;

%debug: visualise grid for the pair
debug = 0;
if debug
    figure(8)
    clf;
    plot(real(rin_pair),imag(rin_pair),'r+')
    hold on
    plot(real(rout_f),imag(rout_f),'bo')
end

project = ~isempty(pair_target_rbm_map);

S = stokSLPmat(rin_pair,rout_f,mu);

if project
    Atot = S + (pair_target_rbm_map - (S*pair_rbm_map)/pair_moment_gram) * ...
        pair_moment_map;
else
    Atot = S;
end

tol = 1e-11; %SVD truncation level
visualise = 0;
[Y,U] = getPseudoFactors(Atot,tol,visualise,svd_opts);

end


