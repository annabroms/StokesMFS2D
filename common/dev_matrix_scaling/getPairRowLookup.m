function pair_lookup = getPairRowLookup(pairs,P)
%GETPAIRROWLOOKUP Build a symmetric lookup from particle pair to pair row.
%
% Syntax:
%   pair_lookup = getPairRowLookup(pairs,P)
%
% Inputs:
%   pairs - N-by-2 close-pair list with one row per unordered pair
%   P     - number of particles
%
% Output:
%   pair_lookup - P-by-P integer matrix. pair_lookup(i,j) is the row index
%                 in pairs for the unordered pair {i,j}, or zero if absent.

if nargin < 2
    error('getPairRowLookup requires pairs and P.');
end

pair_lookup = zeros(P,P);
for row = 1:size(pairs,1)
    i = pairs(row,1);
    j = pairs(row,2);
    pair_lookup(i,j) = row;
    pair_lookup(j,i) = row;
end

end
