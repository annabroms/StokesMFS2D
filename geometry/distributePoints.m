function [z,n,zpp,fz,L] = distributePoints(a,b,tau,N,type)
% DISTRIBUTEPOINTS(a, b, tau, N, type) generates N points evenly spaced in arc length 
% along a closed planar curve. The curve is either an ellipse or a star-shaped curve 
% depending on the input 'type'.
%
% Inputs:
%   a, b   - Shape parameters. For an ellipse, these are the semiaxes (a = x-axis, b = y-axis).
%            For a star, 'a' controls amplitude and 'b' is the frequency.
%   tau     - Rotation
%   N      - Number of points to distribute along the curve.
%   type   - A string specifying the curve type:
%              'ellipse' — standard ellipse with semiaxes a and b.
%              otherwise — star-shaped curve defined by r(t) = (1 + a*sin(b*t)).
%
% Outputs:
%   z      - N array of complex valued point coordinates along the curve.
%   n      - N array of outward unit normal  at each point
%   zpp    - N array of second derivatives at each point
%   fz     - function handle to the generating function for z
%   L      - Total arclength

if nargin<5
    type = 'ellipse';
end

% Parametrize the ellipse
t = chebfun('t', [0, 2*pi]);

if isequal(type,'ellipse')
    x = a * cos(t);
    y = b * sin(t);
else    
    x = (1+a*sin(b*t)).*cos(t);
    y = (1+a*sin(b*t)).*sin(t);   
end

%rotate
z = (x+1i*y)*exp(1i*tau);
x = real(z); 
y = imag(z); 

% Compute arc length as a chebfun
dsdt = sqrt(diff(x).^2 + diff(y).^2);
s_of_t = cumsum(dsdt);             % s(t), arc length as a function of t
L = s_of_t(end);                   % total perimeter of the ellipse

% Invert to get t(s), the inverse mapping
t_of_s = inv(s_of_t);              % now t = t(s)

% Define arc-length parametrized curve
s = linspace(0, L, N);

if isequal(type,'ellipse')
    x_d = @(t) a * cos(t_of_s(t));
    y_d = @(t) b * sin(t_of_s(t));
    fz = @(t) (x_d(t)+1i*y_d(t))*exp(1i*tau);

    z = (x_d(s)+1i*y_d(s))*exp(1i*tau);
    
    ny = -a*sin(t_of_s(s))./sqrt(a^2*sin(t_of_s(s)).^2+b^2*cos(t_of_s(s)).^2);
    nx = -b*cos(t_of_s(s))./sqrt(a^2*sin(t_of_s(s)).^2+b^2*cos(t_of_s(s)).^2);

    n = (nx+1i*ny)*exp(1i*tau);

    %second derivatives
    zpp = -z; %check!
else
    % star shaped domain
    fz = @(t) (1+a*sin(b*t_of_s(t))).*exp(1i*t_of_s(t))*exp(1i*tau);
    fzz = @(t) (a*b*cos(b*t_of_s(t))).*exp(1i*t_of_s(t))*exp(1i*tau)+1i*fz(t);
    fzzz = @(t) ((-a*b^2*sin(b*t_of_s(t)) + 2i*a*b*cos(b*t_of_s(t)) + (1 + a*sin(b*t_of_s(t)))*(-1)) .* exp(1i*t_of_s(t)))*exp(1i*tau);

    z = fz(s).';
    zp = fzz(s).';
    zpp = fzzz(s).';
    n = (-imag(zp)+1i*real(zp))./abs(zp);
end


%debug
% figure()
% subplot(1,2,1)
% plot(x_d,y_d,'.-')
% hold on
% quiver(x_d,y_d,n_x,n_y)
% axis equal

% subplot(1,2,2)
% plot(x,y,'+-')
% plot(a*cos(svals),b*sin(svals),'.-')
% axis equal
% title('Equally spaced points along ellipse')

end


