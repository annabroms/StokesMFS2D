function [xj, w] = gradedptr(N, beta)
    % Function to compute nodes and weights for graded periodic trapezoid rule
    % Supposed to be MATLAB version of Alex's Julia code.
    
    % Define the dx/dt function
    xp = @(t) cosh(beta * sin(t/2));
    
    % Compute derivative values at equispaced nodes on [0, 2pi)
    j = 0:N-1;
    t = 2*pi*j/N;
    xpj = xp(t);
    
    % Compute integral of xp using trapezoid rule (uniform spacing)
    Ix = 2*pi/N * sum(xpj);
    
    % Normalize xpj so total sweep is 2pi
    xpj = xpj * (2*pi / Ix);
    
    % Compute antiderivative using periodic spectral integration
    xj = perispec_antideriv(xpj);
    
    % Return nodes and weights
    w = 2*pi/N * xpj;
end
