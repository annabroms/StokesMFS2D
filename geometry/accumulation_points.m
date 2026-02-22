function [xacc1, xacc2] = accumulation_points(a1, a2, delta, c1, c2)
%ACCUMULATION_POINTS Compute accumulation points for two circular objects.
%
%   [xacc1, xacc2] = accumulation_points(a1, a2, delta, c1, c2)
%
%   INPUTS:
%     a1     - radius of object 1
%     a2     - radius of object 2
%     delta  - distance between objects
%     c1     - 2x1 vector, center of object 1
%     c2     - 2x1 vector, center of object 2
%
%   OUTPUTS:
%     xacc1  - 2x1 vector, accumulation point object 1
%     xacc2  - 2x1 vector, accumulation point object 2

% Ensure inputs are column vectors
c1 = c1(:);
c2 = c2(:);

% Relative unit direction vector
d = (c2 - c1);
d_unit = d / norm(d);

D = a1 + a2 + delta;

% Compute distance offsets
term1 = (D^2 - a2^2 + a1^2) / (2*D);
sqrt1 = sqrt(term1^2 - a1^2);
offset1 = term1 - sqrt1;

term2 = (-D^2 - a2^2 + a1^2) / (2*D);
sqrt2 = sqrt(term2^2 - a2^2);
offset2 = term2 + sqrt2;

% Accumulation points
xacc1 = c1 + offset1 * d_unit;
xacc2 = c2 + offset2 * d_unit;
end
