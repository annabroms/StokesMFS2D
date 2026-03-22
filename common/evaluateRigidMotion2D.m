function u = evaluateRigidMotion2D(q,U,W,rpts)
%EVALUATERIGIDMOTION2D Evaluate planar rigid-body motion on complex points.
%
% Syntax:
%   u = evaluateRigidMotion2D(q,U,W,rpts)
%
% Inputs:
%   q    - Body center (complex scalar).
%   U    - Translational velocity [Ux, Uy].
%   W    - Angular velocity.
%   rpts - Complex target points.
%
% Output:
%   u - Stacked velocity [u_x; u_y] at the target points.

U = U(:);
if numel(U) ~= 2
    error('evaluateRigidMotion2D:badTranslation', ...
        'U must contain exactly two components.');
end

rpts = rpts(:);
u = [U(1) - W*imag(rpts-q); ...
     U(2) + W*real(rpts-q)];

end
