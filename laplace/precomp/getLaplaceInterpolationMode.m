function mode = getLaplaceInterpolationMode(opt)
%GETLAPLACEINTERPOLATIONMODE Validate opt.use_interpolation.
%
% The supported values are:
%   'none'            build an exact map for every separation group;
%   'full'            interpolate the full canonical Cmap in alpha;
%   'reduced'         interpolate Cref+U*B(alpha)*V' panelwise;
%   'reduced_noconst' interpolate U*B(alpha)*V' panelwise (preferred).

value = getOptField(opt,'use_interpolation','none');
if isstring(value)
    if ~isscalar(value)
        error('getLaplaceInterpolationMode:BadValue', ...
            'opt.use_interpolation must be a scalar string or character vector.');
    end
    value = char(value);
end
if ~ischar(value)
    error('getLaplaceInterpolationMode:BadValue', ...
        ['opt.use_interpolation must be ''none'', ''full'', ''reduced'', ', ...
         'or ''reduced_noconst''.']);
end

mode = lower(strtrim(value));
if ~any(strcmp(mode,{'none','full','reduced','reduced_noconst'}))
    error('getLaplaceInterpolationMode:BadValue', ...
        ['opt.use_interpolation must be ''none'', ''full'', ''reduced'', ', ...
         'or ''reduced_noconst''.']);
end
end
