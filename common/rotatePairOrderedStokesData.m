function data_rot = rotatePairOrderedStokesData(data, n, phase, vec_rot)
%ROTATEPAIRORDEREDSTOKESDATA Rotate pair-ordered Stokes circle data.
%
% Input ordering:
%   vector: [x_1; ...; x_M; y_1; ...; y_M] with each block length n
%   matrix: n-by-(2*M) with columns [x_1 ... x_M y_1 ... y_M]

is_vec = isvector(data);
if is_vec
    data = reshape(data(:),n,[]);
end

data_rot = rotateStokesUniformCircleData(data,phase,vec_rot);

if is_vec
    data_rot = data_rot(:);
end
end
