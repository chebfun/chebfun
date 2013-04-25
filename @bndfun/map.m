function f = map(f,newdom)
%MAP   Maps the domain of a BNDFUN F to NEWDOM via a linear change of variable.
%   G = MAP(F,NEWDOM), where the BNDFUN F has a domain [a,b], returns a
%   BNDFUN G defined on [c,d], where c = NEWDOM(1), d = NEWDOM(2), such that
%       G(x) = F(a*(d-x)/(d-c) + b*(x-c)/(d-c)) for all x in [c,d].

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Assign a new linear map to f, obtained from the FUN.LINEAR() method:
f.mapping = fun.linear(newdom);

% Assign the new endpoints of f to the domain property:
f.domain = newdom;

end