function g = cell2mat(f)
%MAT2CELL  Convert an ARRAY of FUNCHEB objects into a vector-valued FUNCHEB.
%
%   G = CELL2MAT(F) converts the FUNCHEB array F in to a single vector-valued
%   FUNCHEB G. F should be a vector array (i.e., not a matrix).
%
% Example:
%   f = funcheb.constructor(@sin);
%   g = funcheb.constructor(@cos);
%   h = cell2mat([f g]);
%
% See also MAT2CELL.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return an empty result:
if ( isempty(f) || numel(f) == 1 )
    g = f;
    return
end

% Extract the values from each of the FUNCHEB objects:
values = {f.values};

% These should all have the same length...
lengths = cellfun(@(x) size(x, 1), values);
minlength = min(lengths);
maxlength = max(lengths);
% If not, then prolong and get the new values:
if ( minlength ~= maxlength )
    for k = 1:numel(f)
        f(k) = prolong(f(k), maxlength);
    end
    values = {f.values};
end

% Append new data to an empty funcheb:
g = f.make(); % Make an empty funcheb.
g.vscale = ([f.vscale]);
g.ishappy = min([f.ishappy]);
g.epslevel = max([f.epslevel]);
g.values = cell2mat(values);
g.coeffs = [f.coeffs];

end