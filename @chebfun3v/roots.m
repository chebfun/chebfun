function varargout = roots(F, varargin)
%ROOTS   Find common zeros of a CHEBFUN3V object.
%   R = ROOTS(F) finds the common zeros of the three CHEBFUN3 objects F(1),
%   F(2) and F(3) in their domain of definition under the assumption that 
%   the solution set is zero-dimensional. Output is a row vector storing 
%   the x-values, y-values, and z-values of the common roots. This script 
%   is also called by the syntax ROOTS(F, G, H), where F, G and H are 
%   CHEBFUN3 objects.
%
% See also CHEBFUN3/ROOT, and CHEBFUN2/ROOTS.

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
        r = roots(f, g, h);
        xroots = r(:, 1); 
        yroots = r(:, 2);
        zroots = r(:, 3);
end

if ( nargout <= 1 )
    varargout{1} = [xroots, yroots, zroots];
else
    varargout = {xroots, yroots, zroots};
end

end