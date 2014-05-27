function f = cos(f)
%COS Cosine of a chebfun2.
%
% See also COSH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

f.chebfun2 = cos(f.chebfun2);

f.der = (-sin(f.chebfun2))*f.der;

end