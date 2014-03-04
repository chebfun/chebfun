function fx = toValues(disc,f)
%TOVALUES  Convert CHEBFUN to a COLLOC1 discretization.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

x = functionPoints(disc);
fx = f(x);

end
