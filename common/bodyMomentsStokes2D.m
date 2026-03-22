function mom = bodyMomentsStokes2D(rsrc,q,lambda_body)
%BODYMOMENTSSTOKES2D Force/torque moments of a bodywise Stokeslet density.
%
% Syntax:
%   mom = bodyMomentsStokes2D(rsrc,q,lambda_body)
%
% Inputs:
%   rsrc        - Complex source locations for one body.
%   q           - Body center.
%   lambda_body  - Stacked [x; y] source strengths for that body.
%
% Output:
%   mom - [Fx; Fy; T] induced by the sources.

rsrc = rsrc(:);
lambda_body = lambda_body(:);
n = numel(rsrc);
if numel(lambda_body) ~= 2*n
    error('bodyMomentsStokes2D:badDensitySize', ...
        'lambda_body must have length 2*numel(rsrc).');
end

lambda_x = lambda_body(1:n);
lambda_y = lambda_body(n+1:end);
dx = real(rsrc - q);
dy = imag(rsrc - q);

mom = [sum(lambda_x); ...
       sum(lambda_y); ...
       sum((-dy).*lambda_x + dx.*lambda_y)];

end
