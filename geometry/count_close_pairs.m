function [n_pairs, pairs, gaps] = count_close_pairs(q, delta_pair, rad)
%COUNT_CLOSE_PAIRS Count particle pairs whose surface gap is below delta_pair.
%
% Syntax:
%   n_pairs = count_close_pairs(q, delta_pair)
%   n_pairs = count_close_pairs(q, delta_pair, rad)
%   [n_pairs, pairs, gaps] = count_close_pairs(...)
%
% Inputs:
%   q          - Complex-valued particle centers (P x 1).
%   delta_pair - Surface-gap threshold.
%   rad        - Optional particle radius/radii. Default: 1.
%
% Outputs:
%   n_pairs - Number of pairs (i,j), i<j, with gap < delta_pair.
%   pairs   - N-by-2 array of close pair indices.
%   gaps    - Gaps corresponding to rows of pairs.
%
% Notes:
%   For particles i and j with radii r_i and r_j, the gap is
%       abs(q(j)-q(i)) - (r_i + r_j).
%
%   Apr 2026

q = q(:);
P = numel(q);

validateattributes(delta_pair,{'numeric'},{'scalar','finite'}, ...
    mfilename,'delta_pair',2);

if nargin < 3 || isempty(rad)
    rad = ones(P,1);
elseif isscalar(rad)
    rad = repmat(rad,P,1);
else
    rad = rad(:);
    if numel(rad) ~= P
        error('count_close_pairs:badRadii', ...
            'rad must be scalar or have one entry per particle.');
    end
end

max_pairs = P*(P-1)/2;
pairs = zeros(max_pairs,2);
gaps = zeros(max_pairs,1);
n_pairs = 0;

for i = 1:P-1
    for j = i+1:P
        gap = abs(q(j)-q(i)) - (rad(i) + rad(j));
        if gap < delta_pair
            n_pairs = n_pairs + 1;
            pairs(n_pairs,:) = [i, j];
            gaps(n_pairs) = gap;
        end
    end
end

pairs = pairs(1:n_pairs,:);
gaps = gaps(1:n_pairs);
end
