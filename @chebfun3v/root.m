function varargout = root(F, varargin)
%ROOTS   Find one common zero of a CHEBFUN3V object.
%   R = ROOT(F) finds one common zero of the three CHEBFUN3 objects F(1), 
%   F(2) and F(3) in their domain of definition under the assumption 
%   that the solution set is zero-dimensional. R is a matrix with three 
%   columns storing the x-value, y-value, and z-value of the common root.
%   This script is also called by the syntax ROOT(F, G, H), where F, G and
%   H are CHEBFUN3 objects.
%
%   [x, y, z] = ROOTS(F) returns the x-value, y-value and z-value as 
%   three separate columns.
%
% See also CHEBFUN3/ROOTS and CHEBFUN2/ROOTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) )
    varargout = {[]};
    return
end

f = F.components{1};
g = F.components{2};
h = F.components{3};

if ( isempty(varargin) )
        r = root(f, g, h);
        xroot = r(:, 1); 
        yroot = r(:, 2);
        zroot = r(:, 3);
end

if ( nargout <= 1 )
    varargout{1} = [xroot, yroot, zroot];
else
    varargout = {xroot, yroot, zroot};
end

end