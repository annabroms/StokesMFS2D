function pair = getLaplacePairInstance(pair_cache,row)
%GETLAPLACEPAIRINSTANCE Recover one actual pair from a cached canonical group.
%
% Syntax:
%   pair = getLaplacePairInstance(pair_cache,row)
%
% Inputs:
%   pair_cache - cache returned by getPairBasisLaplace
%   row       - row index into the actual close-pair list
%
% Output:
%   pair      - struct with fields:
%               meta     : per-pair geometry/rotation metadata
%               group    : shared canonical cached operators for this separation
%               rimage_i : actual image points for the first particle
%               rimage_j : actual image points for the second particle
%
% The cache stores geometry in a canonical frame where the pair midpoint is
% at the origin and the center-to-center line is aligned with the real axis.
% This helper reconstructs the actual pair geometry by applying the inverse
% rigid motion encoded in pair_cache.meta(row).

meta = pair_cache.meta(row);
group = pair_cache.groups(meta.group_id);

pair = struct();
pair.meta = meta;
pair.group = group;
pair.rimage_i = map_points_from_canonical(group.rimage_canon{1},meta);
pair.rimage_j = map_points_from_canonical(group.rimage_canon{2},meta);

end

function z_actual = map_points_from_canonical(z_canon,meta)
%MAP_POINTS_FROM_CANONICAL Map canonical pair-local points back to the solve geometry.

z_canon = z_canon(:);
if isempty(z_canon)
    z_actual = zeros(0,1);
    return
end

z_actual = meta.mid + meta.rot*z_canon;

end
