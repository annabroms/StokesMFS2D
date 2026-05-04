function z = createPeanut(q1, q2, Np, debug, R)
%CREATEPEANUT Creates discretisation of separation surface for two circular
% particles.
%
% Syntax:
%   z = createPeanut(q1, q2, Np, debug)
%   z = createPeanut(q1, q2, Np, debug, R)
%
% Input:
%   q1    - Complex center of circle 1
%   q2    - Complex center of circle 2
%   Np    - Number of discrete points on peanut
%   debug - Boolean to draw peanut
%   R     - Optional circle radius (default 1)

    if nargin < 5 || isempty(R)
        R = 1;
    end

    % --- Geometry ---
    dz    = q2 - q1;
    d     = abs(dz);
    delta = d - 2*R;

    % Rotation angle via angle() — avoids manual acos + sign correction
    theta = angle(dz);

    % Peanut opening half-angle
    alpha = acos((2*R + delta) / (4*R));

    % --- Segment point counts ---
    n_long  = 3*Np/8;       % points on long arcs (t1, t2)
    n_short = Np/8;         % interior points on bridge arcs (t3, t4)

    % --- Parametric angles ---
    t1 = linspace(alpha,        2*pi - alpha,   n_long);
    t2 = linspace(pi + alpha,   3*pi - alpha,   n_long);

    % Endpoint-exclusive linspace: avoids the +2 / trim pattern
    dt3 = (pi  - alpha - alpha) / (n_short + 1);
    dt4 = (pi  - alpha - alpha) / (n_short + 1);
    t3  = alpha + dt3 * (1:n_short);
    t4  = (pi + alpha) + dt4 * (1:n_short);

    % --- Complex arc points using exp(1i*t) — single trig call per segment ---
    e1 = exp(1i * t1);
    e2 = exp(1i * t2);
    e3 = exp(1i * t3);
    e4 = exp(1i * t4);

    % Precompute repeated scalar
    sa = sin(alpha);   % used in z3, z4 imaginary offsets

    z1 = R * e1;
    z2 = (2*R + delta) + R * e2;
    z3 = (R + delta/2) + R*real(e3) + 1i*(-2*R*sa + R*imag(e3));
    z4 = (R + delta/2) + R*real(e4) + 1i*( 2*R*sa + R*imag(e4));

    % --- Concatenate all segments ---
    z = [z1, z2, z3, z4];

    % --- Rotate using complex multiplication (no rotation matrix needed) ---
    rot = exp(1i * theta);   % unit complex rotation
    z   = rot * z;

    % --- Translate to midpoint of q1, q2 ---
    qmid = 0.5 * (q1 + q2);
    z    = z - mean(z) + qmid;

    % Ensure column vector output
    z = z(:);

    % --- Debug plot ---
    if debug
        figure(101)
        subplot(1,2,1)
        plot(real(z1), imag(z1), 'r.--');  hold on
        plot(real(z2), imag(z2), 'b+--');
        plot(real(z3), imag(z3), 'mo--');
        plot(real(z4), imag(z4), 'cs--');
        axis equal;  hold off

        figure(2)
        magenta = [0.8, 0.0, 0.8];
        plot(real(z), imag(z), '.', 'Color', magenta)
        axis equal
    end

end