function phase = getUniformCircleRotationPhase(n,rot)
%GETUNIFORMCIRCLEROTATIONPHASE Fourier phase vector for uniform-circle shifts.
%
% Syntax:
%   phase = getUniformCircleRotationPhase(n,rot)
%
% Inputs:
%   n   - Number of uniform circle samples.
%   rot - Unit complex number exp(1i*theta).
%
% Output:
%   phase - Length-n vector such that ifft(fft(values).*phase) rotates the
%           samples by rot.

if nargin < 1 || isempty(n)
    phase = zeros(0,1);
    return
end

if n == 0
    phase = zeros(0,1);
    return
end

if nargin < 2 || isempty(rot)
    phase = ones(n,1);
    return
end

theta = atan2(imag(rot),real(rot));
if abs(theta) < 10*eps
    phase = ones(n,1);
    return
end

persistent cached_n cached_modes
if isempty(cached_n)
    cached_n = zeros(0,1);
    cached_modes = cell(0,1);
end

idx = find(cached_n == n,1);
if isempty(idx)
    if mod(n,2) == 0
        modes = [0:n/2 -n/2+1:-1];
    else
        modes = [0:(n-1)/2 -(n-1)/2:-1];
    end
    cached_n(end+1,1) = n;
    cached_modes{end+1,1} = modes(:);
    idx = numel(cached_modes);
end

phase = exp(1i*theta*cached_modes{idx});

end
