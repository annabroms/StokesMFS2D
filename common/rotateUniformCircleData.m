function values_rot = rotateUniformCircleData(values,rot)
%ROTATEUNIFORMCIRCLEDATA Periodically shift uniform circle samples by a rigid rotation.
%
% Syntax:
%   values_rot = rotateUniformCircleData(values,rot)
%
% Inputs:
%   values - Samples of a 2*pi-periodic quantity at uniform angles
%            t_k = 2*pi*k/N, stored as N x nd.
%   rot    - Unit complex number exp(1i*theta). The output satisfies
%            values_rot(t_k) = values(t_k + theta).
%
% Output:
%   values_rot - Shifted samples, same size as values.
%
% Notes:
%   This is used to move scalar densities between a fixed proxy-circle
%   parameterization and a pair-local canonical frame.

if nargin < 2 || isempty(rot) || isempty(values)
    values_rot = values;
    return
end

theta = angle(rot);
if abs(theta) < 10*eps
    values_rot = values;
    return
end

n = size(values,1);
if n == 0
    values_rot = values;
    return
end

modes = get_fourier_modes(n);
phase = exp(1i*theta*modes(:));

F = fft(values,[],1);
values_rot = ifft(F.*phase,[],1);

if isreal(values)
    values_rot = real(values_rot);
end

end

function modes = get_fourier_modes(n)
% Return Fourier mode numbers matching MATLAB's fft ordering.

if mod(n,2) == 0
    modes = [0:n/2 -n/2+1:-1];
else
    modes = [0:(n-1)/2 -(n-1)/2:-1];
end

modes = modes(:);

end
