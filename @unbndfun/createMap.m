function m = createMap(ends)
%CREATEMAP   Create a map structure for UNBNDFUN objects.
%   M = CREATEMAP(ENDS) returns a structure that defines a nonlinear map from
%   [-1 1] to the unbounded domain [ENDS(1) ENDS(2)]. The structure MAP consists
%   of three function handles, one string.
%    M.FOR is a function that maps [-1,1] to [ENDS(1) ENDS(2)].
%    M.INV is the inverse map.
%    M.DER is the derivative of the map defined in MAP.FOR.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

m = mapping.unbounded(ends);

end
