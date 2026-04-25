function values_rot = rotateUniformCircleData(values, rot, phase)
%ROTATEUNIFORMCIRCLEDATA Periodically shift uniform circle samples by a
% rigid rotation.
%
% Calling with no arguments runs a self-test.
%
% Syntax:
%   values_rot = rotateUniformCircleData(values, rot)
%   values_rot = rotateUniformCircleData(values, rot, phase)
%
% Inputs:
%   values - Samples of a 2*pi-periodic quantity at uniform angles
%            t_k = 2*pi*k/N, stored as (N x nd).
%   rot    - Unit complex number exp(1i*theta). The output satisfies
%            values_rot(t_k) = values(t_k + theta).
%   phase  - Optional precomputed rotation data:
%            * numeric vector from getUniformCircleRotationPhase(N, rot)
%            * struct spec from getUniformCircleRotationSpec(...)
%              spec.mode must already be lowercase.
%
% Output:
%   values_rot - Shifted samples, same size as values.

%% Self-test
if nargin == 0
    selfTestRotateUniformCircleData();
    return
end

%% Fast exits
if isempty(values) || size(values,1) == 0
    values_rot = values;
    return
end

n = size(values,1);

%% Case A: numeric phase vector passed directly — hot path, zero branches
if nargin >= 3 && ~isempty(phase) && ~isstruct(phase)
    % Ensure column without allocating if already correct shape
    if ~iscolumn(phase)
        phase = phase(:);
    end
    % Inline FFT multiply — no function call overhead
    values_rot = ifft(fft(values,[],1).* phase, [], 1);
    if isreal(values)
        values_rot = real(values_rot);
    end
    return
end

%% Case B: struct spec
if nargin >= 3 && ~isempty(phase) && isstruct(phase)
    spec = phase;
    switch spec.mode          % spec.mode must already be lowercase
        case 'shift'
            shift_steps = 0;
            if isfield(spec,'shift_steps') && ~isempty(spec.shift_steps)
                shift_steps = spec.shift_steps;
            end
            values_rot = circshift(values, -shift_steps, 1);

        case 'fft'
            if numel(spec.phase) ~= n
                error('rotateUniformCircleData:BadSpecSize',...
                    'fft spec phase must have %d entries, got %d.',...
                    n, numel(spec.phase));
            end
            ph = spec.phase(:);
            values_rot = ifft(fft(values,[],1).* ph, [], 1);
            if isreal(values)
                values_rot = real(values_rot);
            end

        case 'oversampled_fft'
            M = spec.n_oversampled;
            if spec.n ~= n
                error('rotateUniformCircleData:BadSpecSize',...
                    'oversampled_fft spec expects %d samples, got %d.',...
                    spec.n, n);
            end
            if M < n || mod(M,n) ~= 0
                error('rotateUniformCircleData:BadOversampledSize',...
                    'n_oversampled must be a positive integer multiple of n.');
            end
            if numel(spec.phase) ~= M
                error('rotateUniformCircleData:BadSpecSize',...
                    'oversampled_fft spec phase must have %d entries, got %d.',...
                    M, numel(spec.phase));
            end
            % Zero-pad spectrum directly — saves one FFT vs interpft
            nd  = size(values,2);
            F   = fft(values,[],1);
            Fup = zeros(M, nd, 'like', F);
            half = floor(n/2);
            Fup(1:half+1,:)       = F(1:half+1,:);
            Fup(M-n+half+2:M,:)   = F(half+2:n,:);
            Fup = Fup * (M/n);
            ph  = spec.phase(:);
            values_up  = ifft(Fup.* ph, [], 1);
            stride     = M/n;
            values_rot = values_up(1:stride:end,:);
            if isreal(values)
                values_rot = real(values_rot);
            end

        otherwise
            error('rotateUniformCircleData:BadMode',...
                'Unsupported rotation mode "%s".', spec.mode);
    end
    return
end

%% Case C: struct passed as rot (legacy convention)
if nargin >= 2 && ~isempty(rot) && isstruct(rot)
    % Recurse once with arguments swapped — only happens at setup, not hot path
    values_rot = rotateUniformCircleData(values, [], rot);
    return
end

%% Case D: plain scalar rot — compute phase then apply inline
if nargin < 2 || isempty(rot)
    values_rot = values;
    return
end

theta = angle(rot);
if abs(theta) < 10*eps
    values_rot = values;
    return
end

phase = getUniformCircleRotationPhase(n, rot);
values_rot = ifft(fft(values,[],1).* phase, [], 1);
if isreal(values)
    values_rot = real(values_rot);
end

end % rotateUniformCircleData

% =========================================================================
%  SELF-TEST  (not on hot path — helper decomposition is fine here)
% =========================================================================

function selfTestRotateUniformCircleData()
fprintf('rotateUniformCircleData self-test\n');
tol      = 1e-10;
all_pass = true;
N        = 64;
nd       = 3;
t        = (0:N-1)' * (2*pi/N);

% Band-limited test signal: exact under all rotation modes
modes  = [1 3 5];
amps   = [1.0 0.5 0.25];
values = zeros(N, nd);
for k = 1:numel(modes)
    values = values + amps(k) * cos(modes(k)*t + (0:nd-1)*0.3);
end

thetas = [pi/7, -pi/3, 2*pi/5, 0, pi];

all_pass = runTest(all_pass, '[1] plain rot', @() testPlainRot(values,thetas,N,tol));
all_pass = runTest(all_pass, '[2] numeric phase', @() testNumericPhase(values,thetas,N,tol));
all_pass = runTest(all_pass, '[3] fft spec', @() testFftSpec(values,thetas,N,tol));
all_pass = runTest(all_pass, '[4] shift spec', @() testShiftSpec(values,N,tol));
all_pass = runTest(all_pass, '[5] oversampled_fft spec', @() testOversampledSpec(values,thetas,N,tol));
all_pass = runTest(all_pass, '[6] identity shortcuts', @() testIdentity(values));
all_pass = runTest(all_pass, '[7] real input stays real', @() testRealPreservation(values));
all_pass = runTest(all_pass, '[8] empty/zero-row inputs', @() testEmptyInputs());

if all_pass
    fprintf('rotateUniformCircleData: all tests passed.\n');
else
    warning('rotateUniformCircleData:SelfTestFailed',...
        'One or more self-test cases failed.');
end
end

function all_pass = runTest(all_pass, label, fn)
fprintf('  %s... ', label);
try
    msg = fn();
    if isempty(msg)
        fprintf('pass\n');
    else
        fprintf('FAIL: %s\n', msg);
        all_pass = false;
    end
catch ex
    fprintf('ERROR: %s\n', ex.message);
    all_pass = false;
end
end

function msg = testPlainRot(values, thetas, N, tol)
msg = '';
for theta = thetas
    rot  = exp(1i*theta);
    vrot = rotateUniformCircleData(values, rot);
    err  = relerr(vrot, referenceRotation(values, theta, N));
    if err > tol
        msg = sprintf('theta=%.4f err=%.2e', theta, err); return
    end
end
end

function msg = testNumericPhase(values, thetas, N, tol)
msg = '';
for theta = thetas
    rot   = exp(1i*theta);
    phase = getUniformCircleRotationPhase(N, rot);
    vrot  = rotateUniformCircleData(values, [], phase);
    err   = relerr(vrot, referenceRotation(values, theta, N));
    if err > tol
        msg = sprintf('theta=%.4f err=%.2e', theta, err); return
    end
end
end

function msg = testFftSpec(values, thetas, N, tol)
msg = '';
for theta = thetas
    rot  = exp(1i*theta);
    spec = struct('mode','fft',...
        'phase', getUniformCircleRotationPhase(N, rot));
    vrot = rotateUniformCircleData(values, [], spec);
    err  = relerr(vrot, referenceRotation(values, theta, N));
    if err > tol
        msg = sprintf('theta=%.4f err=%.2e', theta, err); return
    end
end
end

function msg = testShiftSpec(values, N, tol)
msg = '';
for s = [0 1 3 N-1]
    theta = s * (2*pi/N);
    spec  = struct('mode','shift','shift_steps',s);
    vrot  = rotateUniformCircleData(values, [], spec);
    err   = relerr(vrot, referenceRotation(values, theta, N));
    if err > tol
        msg = sprintf('s=%d err=%.2e', s, err); return
    end
end
end

function msg = testOversampledSpec(values, thetas, N, tol)
msg = '';
M = 4*N;
for theta = thetas
    rot   = exp(1i*theta);
    phase = getUniformCircleRotationPhase(M, rot);
    spec  = struct('mode','oversampled_fft','n',N,'n_oversampled',M,...
        'phase',phase);
    vrot  = rotateUniformCircleData(values, [], spec);
    err   = relerr(vrot, referenceRotation(values, theta, N));
    if err > tol
        msg = sprintf('theta=%.4f err=%.2e', theta, err); return
    end
end
end

function msg = testIdentity(values)
msg = '';
v1 = rotateUniformCircleData(values, 1+0i);
v2 = rotateUniformCircleData(values, exp(1i*0));
v3 = rotateUniformCircleData(values, [], []);
err = max([norm(v1(:)-values(:),inf), norm(v2(:)-values(:),inf),...
           norm(v3(:)-values(:),inf)]);
if err ~= 0
    msg = sprintf('err=%.2e', err);
end
end

function msg = testRealPreservation(values)
msg = '';
vrot = rotateUniformCircleData(real(values), exp(1i*pi/5));
if ~isreal(vrot)
    msg = 'output is complex for real input';
end
end

function msg = testEmptyInputs()
msg = '';
v1 = rotateUniformCircleData([], exp(1i*pi/3));
v2 = rotateUniformCircleData(zeros(0,2), exp(1i*pi/3));
if ~isempty(v1) || ~isempty(v2)
    msg = 'non-empty output for empty input';
end
end

function vref = referenceRotation(values, theta, N)
F     = fft(values, [], 1);
k     = [0:N/2, -N/2+1:-1]';
phase = exp(1i*k*theta);
vref  = ifft(F.* phase, [], 1);
if isreal(values)
    vref = real(vref);
end
end

function e = relerr(a, b)
e = norm(a(:)-b(:), inf) / max(1, norm(b(:), inf));
end