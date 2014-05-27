function f = sin(f) 
% SIN   Sine of a chebfun2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

f.der = cos(f.chebfun2)*f.der;

f.chebfun2 = sin(f.chebfun2);

end