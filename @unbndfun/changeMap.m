function f = changeMap(f, newDom)
%MAP   Map the domain of F to NEWDOM via a nonlinear change of variable.
%   F = CHANGEMAP(F, NEWDOM) returns the original F but defined on a new domain
%   [NEWDOM(1) NEWDOM(2)]. The new domain must remain unbounded.
%
% See also CREATEMAP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Assign a new nonlinear map to f:
f.mapping = unbndfun.createMap(newDom);

% Assign the new endpoints of f to the domain property:
f.domain = newDom;

end
