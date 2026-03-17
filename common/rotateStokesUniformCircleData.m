function data_rot = rotateStokesUniformCircleData(data, phase, vec_rot)
%ROTATESTOKESUNIFORMCIRCLEDATA Rotate vector-valued circle data.
%
% Syntax:
%   data_rot = rotateStokesUniformCircleData(data, phase, vec_rot)
%
% Input layout:
%   data is N-by-(2*M) with columns [x_1 ... x_M y_1 ... y_M].
%
% The rotation consists of:
%   1) a sample shift on each uniform circle via rotateUniformCircleData
%   2) a 2D vector rotation by vec_rot (complex unit number)

if nargin < 2
    phase = [];
end
if nargin < 3 || isempty(vec_rot)
    vec_rot = 1;
end

if isempty(data)
    data_rot = data;
    return
end

[n, m2] = size(data);
if mod(m2,2) ~= 0
    error('rotateStokesUniformCircleData expects [x ... y ...] columns.');
end

m = m2/2;
z = data(:,1:m) + 1i*data(:,m+1:end);
z = rotateUniformCircleData(z,[],phase);
z = vec_rot*z;

data_rot = zeros(n,m2);
data_rot(:,1:m) = real(z);
data_rot(:,m+1:end) = imag(z);
end
