function L = bc(L, c)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
%BC  Set linop constraints (overwrite existing).
validateattributes(c, {'linopConstraint'})
L.constraint = c;
end
