function F = perispec_antideriv(f)
    % perispec_antideriv computes periodic spectral antiderivative of f
    % f should be a complex (or real) vector
    % returns complex antiderivative F with mean(f) * linear ramp added
    % Supposed to be MATLAB version of Alex's Julia code.
    
    N = length(f);
    
    % Construct the spectral filter
    n = [0, 1:floor((N-1)/2), -floor(N/2):-1];
    filter = zeros(1, N);
    filter(2:end) = 1 ./ (1i * n(2:end));  % avoid division by zero at n=0
    
    % Compute FFT-based antiderivative (ignores mean)
    F = ifft(fft(f) .* filter);
    
    % Compute and add the sawtooth term to account for the mean
    fbar = mean(f);
    t = (0:N-1) * (2*pi/N);
    F = F + fbar * t - F(1);  % fix offset so F starts at 0
end
