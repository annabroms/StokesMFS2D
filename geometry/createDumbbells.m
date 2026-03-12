function qnew = createDumbbells(P,delta,R,debug_plot)
%CREATEDUMBBELLS Build a geometry of paired disks ("dumbbells").
%
% Syntax:
%   q = createDumbbells(P,delta)
%   q = createDumbbells(P,delta,R)
%   q = createDumbbells(P,delta,R,debug_plot)
%
% Inputs:
%   P          - Number of disks (must be even).
%   delta      - Gap parameter within each dumbbell pair (relative to R).
%   R          - Disk radius. Default: 1.
%   debug_plot - Logical flag. If true, visualise:
%                1) disks of radius R centered at q,
%                2) enclosing pair circles of radius 2R for each dumbbell.
%                Default: false.
%
% Output:
%   qnew       - Complex centers (P x 1).
%
% See also: grow_cluster.
%
% Anna Broms, Mar 2026

if nargin<1
    self_test()
    return; 
end

if nargin < 3 || isempty(R)
    R = 1;
end
if nargin < 4 || isempty(debug_plot)
    debug_plot = false;
end

if mod(P,2) ~= 0
    error('P must be even')
end

% Create centers for P/2 dumbbells, then place each pair with random orientation.
[q,~] = grow_cluster(P/2,R*delta/(2*R+delta*R/2),2,2*R+delta*R/2);
%[q,~] = grow_cluster(P/2,R*delta/(2*R),2,2*R);
qnew = zeros(P,1);
disp('Dumbbell positions generated')

for k = 1:P/2
    theta = 2*pi*rand(1);
    qtemp = (2*R+delta*R)*exp(1i*theta);
    q1 = -qtemp/2;
    q2 = qtemp/2;

    qnew(2*(k-1)+1) = q(k)+q1;
    qnew(2*k) = q(k)+q2;
end

if debug_plot
    t = linspace(0,2*pi,300);
    figure('Name','createDumbbells debug','Color','w');
    clf;
    hold on;

    % Plot all body circles of radius R.
    for i = 1:P
        x = real(qnew(i)) + R*cos(t);
        y = imag(qnew(i)) + R*sin(t);
        plot(x,y,'k-','LineWidth',1.2);
        plot(real(qnew(i)),imag(qnew(i)),'k.','MarkerSize',10);
    end

    % Plot one enclosing circle of radius 2R for each pair.
    for k = 1:P/2
        cpair = q(k);
        xpair = real(cpair) + (2*R+R*delta)*cos(t);
        ypair = imag(cpair) + (2*R+R*delta)*sin(t);
        plot(xpair,ypair,'r--','LineWidth',1.2);
    end

    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');
    title(sprintf('createDumbbells debug (P=%d, delta=%.3g, R=%.3g)',P,delta,R), ...
        'Interpreter','none');
    hold off;
    drawnow;
end

end

function self_test()
% Draw a configuration
    P = 4;
    delta = 0.5;
    R = 2;
    q = createDumbbells(P,delta,R,1);
end

