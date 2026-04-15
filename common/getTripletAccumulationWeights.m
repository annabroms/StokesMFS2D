function weights = getTripletAccumulationWeights(triplets,P,mode)
%GETTRIPLETACCUMULATIONWEIGHTS Construct bodywise weights for triplet accumulation.
%
% Syntax:
%   weights = getTripletAccumulationWeights(triplets,P)
%   weights = getTripletAccumulationWeights(triplets,P,mode)
%
% Inputs:
%   triplets - N-by-3 unordered triplets [i j k]
%   P        - number of particles
%   mode     - 'sum' | 'edge_weighted' | 'pair_subtracted'
%
% Output:
%   weights  - N-by-3 bodywise accumulation weights

if nargin < 2
    error('getTripletAccumulationWeights requires triplets and P.');
end

if nargin < 3 || isempty(mode)
    mode = 'edge_weighted';
end

n_trip = size(triplets,1);
if n_trip == 0
    weights = zeros(0,3);
    return
end

mode = lower(char(mode));
switch mode
    case 'sum'
        weights = ones(n_trip,3);

    case 'edge_weighted'
        edge_mult = zeros(P,P);
        for row = 1:n_trip
            ids = sort(triplets(row,:));
            edge_mult(ids(1),ids(2)) = edge_mult(ids(1),ids(2)) + 1;
            edge_mult(ids(1),ids(3)) = edge_mult(ids(1),ids(3)) + 1;
            edge_mult(ids(2),ids(3)) = edge_mult(ids(2),ids(3)) + 1;
        end

        weights = zeros(n_trip,3);
        for row = 1:n_trip
            ids = sort(triplets(row,:));
            mij = edge_mult(ids(1),ids(2));
            mik = edge_mult(ids(1),ids(3));
            mjk = edge_mult(ids(2),ids(3));

            weights(row,1) = 0.5*(1/max(mij,1) + 1/max(mik,1));
            weights(row,2) = 0.5*(1/max(mij,1) + 1/max(mjk,1));
            weights(row,3) = 0.5*(1/max(mik,1) + 1/max(mjk,1));
        end

    case 'pair_subtracted'
        % Weights are ignored in the hybrid 2B + 3B-remainder mode, but
        % return ones so accidental callers do not fail.
        weights = ones(n_trip,3);

    otherwise
        error('Unsupported triplet_accumulation_mode "%s".',mode);
end

end
