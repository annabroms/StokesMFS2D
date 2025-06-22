function [z, n, zp, zpp] = getPoints(a, b, tau, N, type)
%GETPOINTS Parametrize and discretize a 2D closed curve.
%
%   [z, n, zp, zpp] = GETPOINTS(a, b, tau, N, type) returns N points z on a
%   closed planar curve, together with the complex unit normal vector n,
%   the first derivative zp, and the second derivative zpp with respect to
%   the parametrization.
%
%   INPUTS:
%     a     - Shape parameter: semi-axis for ellipse, amplitude for star.
%     b     - Shape parameter: semi-axis for ellipse, frequency for star.
%     tau   - Rotation angle in radians.
%     N     - Number of discretization points.
%     type  - (Optional) 'ellipse' (default) or 'star' for star-shaped domain.
%
%   OUTPUTS:
%     z     - Complex coordinates of curve points (N×1).
%     n     - Complex unit normal vectors (N×1), outward-pointing.
%     zp    - First derivatives dz/dt (N×1).
%     zpp   - Second derivatives d²z/dt² (N×1).
%
%   NOTES:
%     - For 'star', the curve is given by (1 + a sin(bt)) * exp(it).
%


if nargin<5
    type = 'ellipse';
end

s = linspace(0, 2*pi, N+1);
s = s(1:end-1);

if isequal(type,'ellipse')
    x_d = @(t) a * cos(t);
    y_d = @(t) b * sin(t);   

    z = (x_d(s)+1i*y_d(s))*exp(1i*tau);

    zpx = -a*sin(s);
    zpy = b*cos(s);
    zp = zpx+1i*zpy;
    
    ny = -a*sin(s)./sqrt(a^2*sin(s).^2+b^2*cos(s).^2);
    nx = -b*cos(s)./sqrt(a^2*sin(s).^2+b^2*cos(s).^2);

    n = (nx+1i*ny)*exp(1i*tau);

    %second derivatives
    zpp = -z; %check!
else
    % star shaped domain
    fz = @(t) (1+a*sin(b*t)).*exp(1i*t)*exp(1i*tau);
    fzz = @(t) (a*b*cos(b*t)).*exp(1i*t)*exp(1i*tau)+1i*fz(t);
    fzzz = @(t) ((-a*b^2*sin(b*t) + 2i*a*b*cos(b*t) + (1 + a*sin(b*t))*(-1)) .* exp(1i*t))*exp(1i*tau);

    z = fz(s).';
    zp = fzz(s).';
    zpp = fzzz(s).';
    n = (-imag(zp)+1i*real(zp))./abs(zp);
end


end


