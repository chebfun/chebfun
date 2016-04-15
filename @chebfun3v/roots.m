function varargout = roots( F, varargin )
%ROOTS   Find the common zeros of a CHEBFUN3V object.
%   r = ROOTS(F) finds the common zeros of the three trivariate functions 
%   F(1), F(2) and F(3) in their domain of definition under the assumption 
%   that the solution set is zero-dimensional. R is a matrix with three 
%   columns storing the x-values, y-values, and z-values of the solutions. 
%   This script is also called by the syntax ROOTS(f,g,h), where f, g and
%   h are CHEBFUN3 objects.
%
%   [x, y, z] = ROOTS(F) returns the x-values, y-values and z-values as 
%   three separate columns.
%
% See also CHEBFUN3/ROOTS, and CHEBFUN3/ROOTS.

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
        [xroots, yroots, zroots] = chebfun3.roots(f,g,h);
        xroots = xroots.'; 
        yroots = yroots.';
        zroots = zroots.';
end

if ( nargout <= 1 )
    varargout{1} = [xroots ; yroots, zroots].';
else
    varargout = {xroots, yroots, zroots};
end

end