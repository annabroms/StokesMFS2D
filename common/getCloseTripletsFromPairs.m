function triplets = getCloseTripletsFromPairs(pairs,P)
%GETCLOSETRIPLETSFROMPAIRS Enumerate unordered close triplets from a pair list.
%
% Syntax:
%   triplets = getCloseTripletsFromPairs(pairs,P)
%
% Inputs:
%   pairs - N-by-2 array of close pairs [i j] with i < j
%   P     - number of particles
%
% Output:
%   triplets - M-by-3 array of unordered close triplets [i j k] with
%              i < j < k. A triplet is included iff all three edges are
%              present in the close-pair graph.

if nargin < 2
    error('getCloseTripletsFromPairs requires pairs and P.');
end

if isempty(pairs)
    triplets = zeros(0,3);
    return
end

adj = false(P,P);
for row = 1:size(pairs,1)
    i = pairs(row,1);
    j = pairs(row,2);
    adj(i,j) = true;
    adj(j,i) = true;
end

triplets = zeros(0,3);
for i = 1:P-2
    neigh_i = find(adj(i,i+1:end)) + i;
    for jj = 1:numel(neigh_i)
        j = neigh_i(jj);
        neigh_j = find(adj(j,j+1:end)) + j;
        common = intersect(neigh_i(neigh_i > j),neigh_j,'stable');
        if ~isempty(common)
            triplets = [triplets; [i*ones(numel(common),1), j*ones(numel(common),1), common(:)]];
        end
    end
end

end
