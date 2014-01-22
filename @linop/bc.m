function L = bc(L, c)
%BC  Set linop constraints (overwrite existing).

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

validateattributes(c, {'linopConstraint'})
L.constraint = c;
end
