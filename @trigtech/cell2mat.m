function g = cell2mat(f)
%CELL2MAT   Convert an array of TRIGTECH objects into an array-valued TRIGTECH.
%
%   G = CELL2MAT(F) converts the TRIGTECH array F into a single array-valued
%   TRIGTECH G. F should be a vector array (i.e., not a matrix).
%
% Example:
%   f = trigtech.constructor(@sin);
%   g = trigtech.constructor(@cos);
%   h = cell2mat([f g]);
%
% See also MAT2CELL.
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: This function is probably not needed anymore.

% Return an empty result:
if ( isempty(f) || numel(f) == 1 )
    g = f;
    return
end

% Extract the values from each of the TRIGTECH objects:
values = { f.values };

% These should all have the same length...
lengths = cellfun(@(x) size(x, 1), values);
minlength = min(lengths);
maxlength = max(lengths);

% If not, then prolong and get the new values:
if ( minlength ~= maxlength )
    for k = 1:numel(f)
        f(k) = prolong(f(k), maxlength);
    end
    values = { f.values };
end

% Append new data to an empty TRIGTECH:
g = f.make(); % Make an empty TRIGTECH.
g.ishappy = min([f.ishappy]);
g.values = cell2mat(values);
g.coeffs = [ f.coeffs ];
g.isReal = [ f.isReal ];

end
