function f = map(f,newdom)
%MAP   Maps the domain of an UNBNDFUN F to NEWDOM via a change of variable.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( ~all(isinf(f.domain) == isinf(newdom)) )
    error('CHEBFUN:map:unbnd', 'Domain must remain unbounded.');
end

% Assign a new nonlinear map to f, obtained from the FUN.UNBOUNDED() method:
f.mapping = fun.unbounded(newdom);

% Assign the new endpoints of f to the domain property:
f.domain = newdom;

end