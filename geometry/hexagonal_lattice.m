function q = hexagonal_lattice(delta,R,rad)
%HEXAGONAL_LATTICE Construct a hexagonal disk packing with R rings around 
% the center disk.
%
% Syntax:
%   q = hexagonal_lattice(delta,R)
%   q = hexagonal_lattice(delta,R,rad)
%
% Inputs:
%   delta - Relative gap so nearest-neighbour spacing is rad*(2+delta).
%   R     - Number of hexagonal rings around the center particle.
%   rad   - Disk radius. Default is 1.
%
% Output:
%   q     - Complex-valued particle centers.
%
% Notes:
%   The total number of particles is 1 + 3*R*(R+1).

if nargin == 0
    test_hexagonal_lattice();
    return
end

if nargin < 3 || isempty(rad)
    rad = 1;
end

R = round(R);
if R < 0
    error('hexagonal_lattice:InvalidRingCount','R must be nonnegative.');
end

d = rad*(2 + delta);   % nearest-neighbour distance
a1 = d*[1; 0];
a2 = d*[0.5; sqrt(3)/2];

q = zeros(1 + 3*R*(R+1),1);
idx = 1;
for m = -R:R
    for n = -R:R
        % Hex distance in axial coordinates is (|m|+|n|+|m+n|)/2.
        if abs(m) + abs(n) + abs(m+n) <= 2*R
            r = m*a1 + n*a2;
            q(idx) = r(1) + 1i*r(2);
            idx = idx + 1;
        end
    end
end

q = q(1:idx-1);

end

function test_hexagonal_lattice
delta = 1e-3;
R = 8;
rad = 2;
q = hexagonal_lattice(delta,R,rad);

t = linspace(0,2*pi,240);
figure('Name','hexagonal_lattice self-test','Color','w');
hold on
for k = 1:numel(q)
    zk = q(k) + rad*(cos(t)+1i*sin(t));
    plot(real(zk),imag(zk),'k-','LineWidth',1.1);
    plot(real(q(k)),imag(q(k)),'k.','MarkerSize',10);
end
axis equal
grid on
xlabel('x')
ylabel('y')
title(sprintf('hexagonal\\_lattice self-test: R = %d, rad = %.1f, P = %d', ...
    R,rad,numel(q)),'Interpreter','none')

end
