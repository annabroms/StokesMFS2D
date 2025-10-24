function q = hexagonal_lattice(delta, R)
    % Construct positions for a hexagonal lattice with R rings around the center.
    % Total number of particles: 1 + 3*R*(R+1)
    %
    % q : complex-valued positions of particles in the plane.

    d = 2 + delta;   % nearest neighbor distance
    a1 = d*[1; 0];
    a2 = d*[0.5; sqrt(3)/2];

    q = [];
    for m = -R:R
        for n = -R:R
            if abs(m) + abs(n) + abs(m+n) <= R
                r = m*a1 + n*a2;
                q = [q; r(1) + 1i*r(2)];
            end
        end
    end
end