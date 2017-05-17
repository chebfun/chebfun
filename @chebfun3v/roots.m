function varargout = roots(F)
%ROOTS   Find common zeros of a CHEBFUN3V object.
%   R = ROOTS(F) finds the common zeros of the three CHEBFUN3 objects F(1),
%   F(2) and F(3) in their domain of definition under the assumption that 
%   the solution set is zero-dimensional. R is a matrix with three columns
%   storing the x-values, y-values, and z-values of the common roots. This
%   function is also called by the syntax ROOTS(F, G, H), where F, G and H
%   are CHEBFUN3 objects.
%
%   [x, y, z] = ROOTS(F) returns the x-values, y-values and z-values as
%   three separate columns.
%
% See also CHEBFUN3/ROOTS, and CHEBFUN2/ROOTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) )
    varargout = {[]};
    return
end

f = F.components{1}; 
g = F.components{2};
h = F.components{3};

r = roots(f, g, h);

if ( nargout <= 1 )
    varargout{1} = r;
else
    varargout = {r(:, 1), r(:, 2), r(:, 3)};
end

end