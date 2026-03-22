function geom_check = buildStokes1BCheckGeom(geom,rcheck_b)
%BUILDSTOKES1BCHECKGEOM Create a 1B geometry struct for check-boundary tests.
%
% Syntax:
%   geom_check = buildStokes1BCheckGeom(geom,rcheck_b)
%
% Inputs:
%   geom     - Geometry struct returned by prepareStokes1BEnhanced.
%   rcheck_b - Stacked check-boundary points, bodywise ordered.
%
% Output:
%   geom_check - Copy of geom with the target grid replaced by rcheck_b.

if nargin < 2
    error('buildStokes1BCheckGeom requires geom and rcheck_b.');
end

rcheck_b = rcheck_b(:);
P = numel(geom.q);
if mod(numel(rcheck_b),P) ~= 0
    error('buildStokes1BCheckGeom:badCheckGrid', ...
        'rcheck_b must contain the same number of points for each body.');
end

n_bound = numel(rcheck_b)/P;

geom_check = geom;
geom_check.rout = rcheck_b;
geom_check.rcheck = rcheck_b;
geom_check.rout_body = cell(P,1);
geom_check.target_count = n_bound*ones(P,1);
geom_check.target_offsets = cumsum([1; geom_check.target_count(1:end-1)]);
geom_check.total_target_count = numel(rcheck_b);
geom_check.target_x_ind = cell(P,1);
geom_check.target_y_ind = cell(P,1);

for k = 1:P
    inds = (k-1)*n_bound+1:k*n_bound;
    geom_check.rout_body{k} = rcheck_b(inds);
    geom_check.target_x_ind{k} = inds;
    geom_check.target_y_ind{k} = geom_check.total_target_count + inds;
end

end
