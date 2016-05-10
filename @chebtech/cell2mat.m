function g = cell2mat(f)
%CELL2MAT   Convert an array of CHEBTECH objects into an array-valued CHEBTECH.
%   G = CELL2MAT(F) converts the CHEBTECH array F into a single array-valued
%   CHEBTECH G. F should be a vector array (i.e., not a matrix).
%
% Example:
%   f = chebtech.constructor(@sin);
%   g = chebtech.constructor(@cos);
%   h = cell2mat([f g]);
%
% See also MAT2CELL.
%
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: This function is probably not needed anymore.

% Return an empty result:
if ( isempty(f) || numel(f) == 1 )
    g = f;
    return
end

% Extract the coeffs from each of the CHEBTECH objects:
coeffs = { f.coeffs }; 

% These should all have the same length...
lengths = cellfun(@(x) size(x, 1), coeffs);
minlength = min(lengths);
maxlength = max(lengths);

% If not, then prolong and get the new values:
if ( minlength ~= maxlength )
    for k = 1:numel(f)
        f(k) = prolong(f(k), maxlength);
    end
    coeffs = { f.coeffs };
end

% Append new data to an empty CHEBTECH:
g = f.make(); % Make an empty CHEBTECH.
g.ishappy = min([f.ishappy]);
g.coeffs = cell2mat(coeffs);

end
