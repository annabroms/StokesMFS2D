function values_rot = rotateUniformCircleData(values,rot,phase)
%ROTATEUNIFORMCIRCLEDATA Periodically shift uniform circle samples by a rigid rotation.
%
% Syntax:
%   values_rot = rotateUniformCircleData(values,rot)
%   values_rot = rotateUniformCircleData(values,rot,phase)
%
% Inputs:
%   values - Samples of a 2*pi-periodic quantity at uniform angles
%            t_k = 2*pi*k/N, stored as N x nd.
%   rot    - Unit complex number exp(1i*theta). The output satisfies
%            values_rot(t_k) = values(t_k + theta).
%   phase  - Optional precomputed Fourier phase vector from
%            getUniformCircleRotationPhase(size(values,1),rot).
%
% Output:
%   values_rot - Shifted samples, same size as values.
%
% Notes:
%   This is used to move scalar densities between a fixed proxy-circle
%   parameterization and a pair-local canonical frame.

if isempty(values)
    values_rot = values;
    return
end

n = size(values,1);
if n == 0
    values_rot = values;
    return
end

if nargin < 3 || isempty(phase)
    if nargin < 2 || isempty(rot)
        values_rot = values;
        return
    end
    theta = atan2(imag(rot),real(rot));
    if abs(theta) < 10*eps
        values_rot = values;
        return
    end
    phase = getUniformCircleRotationPhase(n,rot);
else
    phase = phase(:);
    if numel(phase) ~= n
        error('rotateUniformCircleData:BadPhaseSize', ...
            'phase must have one entry per sample.');
    end
end

F = fft(values,[],1);
values_rot = ifft(F.*phase,[],1);

if isreal(values)
    values_rot = real(values_rot);
end

end
