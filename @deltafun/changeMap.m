function f = changeMap(f, newdom)
%CHANGEMAP   Map the domain of a DELTAFUN via a linear change of variable.
%   G = MAP(F,NEWDOM), where the DELTAFUN F has a domain [a, b], returns a
%   DELTAFUN G defined on [c, d], where c = NEWDOM(1), d = NEWDOM(2), such that
%       G(x) = F(a*(d - x)/(d - c) + b*(x - c)/(d - c)) for all x in [c, d].

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Store the old mapping:
oldMapping = f.funPart.mapping;

% Map the funPart:
f.funPart = changeMap(f.funPart, newdom);

% Grab the new mapping:
newMapping = f.funPart.mapping;

% Map the deltaLocs:
f.deltaLoc = newMapping.For(oldMapping.Inv(f.deltaLoc));

end
