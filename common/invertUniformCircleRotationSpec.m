function spec_inv = invertUniformCircleRotationSpec(spec)
%INVERTUNIFORMCIRCLEROTATIONSPEC Build the inverse of a uniform-circle rotation spec.

if isempty(spec)
    spec_inv = spec;
    return
end

if ~isstruct(spec)
    error('invertUniformCircleRotationSpec:BadInput', ...
        'Expected a rotation spec struct.');
end

spec_inv = spec;
if isfield(spec_inv,'rot') && ~isempty(spec_inv.rot)
    spec_inv.rot = conj(spec_inv.rot);
end

switch lower(char(spec.mode))
    case 'shift'
        if isfield(spec,'shift_steps') && ~isempty(spec.shift_steps)
            spec_inv.shift_steps = -spec.shift_steps;
        end

    case 'fft'
        spec_inv.phase = conj(spec.phase);

    case 'oversampled_fft'
        spec_inv.phase = conj(spec.phase);
        spec_inv.shift_steps = 0;

    otherwise
        error('invertUniformCircleRotationSpec:BadMode', ...
            'Unsupported rotation mode "%s".',spec.mode);
end

end
