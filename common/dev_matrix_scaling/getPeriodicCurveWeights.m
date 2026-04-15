function w = getPeriodicCurveWeights(z,center)
%GETPERIODICCURVEWEIGHTS Approximate local arclength weights on a closed curve.
%
% Syntax:
%   w = getPeriodicCurveWeights(z)
%   w = getPeriodicCurveWeights(z,center)
%
% Inputs:
%   z      - Complex boundary nodes on one closed component.
%            If center is omitted, z is assumed to already be ordered
%            around the curve.
%   center - Optional complex center used to sort circular-component nodes
%            by polar angle before the periodic neighbor-distance weights
%            are formed. Returned weights match the original ordering.
%
% Output:
%   w      - Positive local arclength weights, one per node in z.
%
% Anna Broms, Mar 2026

z = z(:);
n = numel(z);

if n == 0
    w = zeros(0,1);
    return
end

if n == 1
    w = 1;
    return
end

if nargin > 1 && ~isempty(center)
    theta = mod(angle(z-center),2*pi);
    [~,order] = sort(theta);
    z_ordered = z(order);
    w_ordered = compute_periodic_weights(z_ordered);
    w = zeros(n,1);
    w(order) = w_ordered;
else
    w = compute_periodic_weights(z);
end

end

function w = compute_periodic_weights(z)
n = numel(z);
prev = z([n; (1:n-1)']);
next = z([(2:n)'; 1]);
w = 0.5*(abs(z-prev) + abs(next-z));
w = real(w);

if any(~isfinite(w))
    error('getPeriodicCurveWeights:InvalidWeights', ...
        'Computed non-finite periodic curve weights.');
end

positive = w > 0;
if ~all(positive)
    if any(positive)
        floor_w = eps(max(1,max(w(positive))));
    else
        floor_w = 1;
    end
    w(~positive) = floor_w;
end
end
