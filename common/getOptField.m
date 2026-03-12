function value = getOptField(opt, name, default_value)
%GETOPTFIELD Return option field value or default when missing/empty.
%
% Syntax:
%   value = getOptField(opt, name, default_value)
%
% Inputs:
%   opt           - Options struct.
%   name          - Field name to retrieve.
%   default_value - Value used when field is absent or empty.
%
% Output:
%   value - Retrieved option value.
%
% See also: isfield.
%
% Anna Broms, Mar 2026

if isstruct(opt) && isfield(opt,name) && ~isempty(opt.(name))
    value = opt.(name);
else
    value = default_value;
end

end
