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
%   phase  - Optional precomputed rotation data. This can be either:
%            * a Fourier phase vector from
%              getUniformCircleRotationPhase(size(values,1),rot)
%            * a rotation spec from getUniformCircleRotationSpec(...)
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

spec = [];

if nargin < 3 || isempty(phase)
    if nargin < 2 || isempty(rot)
        values_rot = values;
        return
    end
    if isstruct(rot)
        spec = rot;
    else
        theta = atan2(imag(rot),real(rot));
        if abs(theta) < 10*eps
            values_rot = values;
            return
        end
        phase = getUniformCircleRotationPhase(n,rot);
    end
elseif isstruct(phase)
    spec = phase;
else
    phase = phase(:);
    if numel(phase) ~= n
        error('rotateUniformCircleData:BadPhaseSize', ...
            'phase must have one entry per sample.');
    end
end

if ~isempty(spec)
    switch lower(char(spec.mode))
        case 'shift'
            shift_steps = 0;
            if isfield(spec,'shift_steps') && ~isempty(spec.shift_steps)
                shift_steps = spec.shift_steps;
            end
            values_rot = circshift(values,-shift_steps,1);

        case 'fft'
            phase = spec.phase(:);
            if numel(phase) ~= n
                error('rotateUniformCircleData:BadSpecSize', ...
                    'rotation spec phase must have one entry per sample.');
            end

            F = fft(values,[],1);
            values_rot = ifft(F.*phase,[],1);

        case 'oversampled_fft'
            if spec.n ~= n
                error('rotateUniformCircleData:BadSpecSize', ...
                    'rotation spec expects %d samples, got %d.',spec.n,n);
            end
            if spec.n_oversampled < n || mod(spec.n_oversampled,n) ~= 0
                error('rotateUniformCircleData:BadOversampledSize', ...
                    'oversampled grid size must be an integer multiple of n.');
            end

            values_up = interpft(values,spec.n_oversampled,1);
            phase = spec.phase(:);
            if numel(phase) ~= spec.n_oversampled
                error('rotateUniformCircleData:BadSpecSize', ...
                    'oversampled rotation spec phase must have one entry per oversampled sample.');
            end
            values_up = ifft(fft(values_up,[],1).*phase,[],1);
            stride = spec.n_oversampled/n;
            values_rot = values_up(1:stride:end,:);

        otherwise
            error('rotateUniformCircleData:BadMode', ...
                'Unsupported rotation mode "%s".',spec.mode);
    end
else
    F = fft(values,[],1);
    values_rot = ifft(F.*phase,[],1);
end

if isreal(values)
    values_rot = real(values_rot);
end

end
