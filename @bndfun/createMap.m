function map = createMap(ends)
%CREATEMAP   Creates a linear map structure for BNDFUN objects.
%   MAP = CREATEMAP(ENDS), where ENDS is a two-vector, returns a structure that
%   defines a linear map. The structure MAP consists of three function handles:
%      MAP.FOR is a function that maps [-1,1] to [ENDS(1), ENDS(2)].
%      MAP.INV is the inverse map.
%      MAP.DER is the derivative of the map defined in MAP.FOR

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

map = mapping.linear(ends);
return

a = ends(1);
b = ends(2);
map = struct('for', @(y) b*(y + 1)/2 + a*(1 - y)/2, ...
             'inv', @(x) (x - a)/(b - a) - (b - x)/(b - a), ...
             'der', @(y) (b - a)/2 + 0*y);
    
% Note: writing the map in this form ensures that -1 is mapped to a and 1
% is mapped to b (in the presence of rounding errors).

end
