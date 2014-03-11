function fx = toValues(disc,f)
%TOVALUES  Convert CHEBFUN to a COLLOC1 discretization.
%   TOVALUES(DISC,F) converts a chebfun F to values at 1st kind points in
%   the DISC.DOMAIN.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

x = functionPoints(disc);
fx = f(x);

end
