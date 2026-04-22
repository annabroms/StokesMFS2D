function spec = getUniformCircleRotationSpec(n,rot,opt)
%GETUNIFORMCIRCLEROTATIONSPEC Build reusable rotation data for uniform-circle shifts.
%
% Syntax:
%   spec = getUniformCircleRotationSpec(n,rot)
%   spec = getUniformCircleRotationSpec(n,rot,opt)
%
% Inputs:
%   n   - Number of uniform circle samples.
%   rot - Unit complex number exp(1i*theta).
%   opt - Optional struct with fields:
%         rotation_mode       : 'fft' | 'oversampled_fft'
%         rotation_oversample : positive integer factor for oversampled_fft
%                               (used to build an oversampled Fourier shift)
%
% Output:
%   spec - Struct consumed by rotateUniformCircleData.

if nargin < 1 || isempty(n) || n == 0
    spec = struct('mode','fft','n',0,'rot',1,'phase',zeros(0,1), ...
        'n_oversampled',0,'shift_steps',0);
    return
end

if nargin < 2 || isempty(rot)
    rot = 1;
end

mode = 'fft';
oversample = 8;

if nargin >= 3 && ~isempty(opt)
    if isstruct(opt)
        mode = lower(char(getOptField(opt,'rotation_mode','fft')));
        oversample = getOptField(opt,'rotation_oversample',oversample);
    else
        mode = lower(char(opt));
    end
end

if isempty(mode)
    mode = 'fft';
end

theta = atan2(imag(rot),real(rot));
shift_steps = round(theta*n/(2*pi));

spec = struct();
spec.mode = mode;
spec.n = n;
spec.rot = rot;
spec.phase = [];
spec.n_oversampled = n;
spec.shift_steps = shift_steps;

if abs(theta - 2*pi*shift_steps/n) <= 100*eps(max(1,abs(theta)))
    spec.mode = 'shift';
    return
end

switch mode
    case 'fft'
        spec.phase = getUniformCircleRotationPhase(n,rot);

    case 'oversampled_fft'
        oversample = round(oversample);
        if oversample < 1
            error('getUniformCircleRotationSpec:BadOversample', ...
                'rotation_oversample must be a positive integer.');
        end
        spec.n_oversampled = oversample*n;
        spec.phase = getUniformCircleRotationPhase(spec.n_oversampled,rot);

    otherwise
        error('getUniformCircleRotationSpec:BadMode', ...
            'Unsupported rotation_mode "%s".',mode);
end

end
