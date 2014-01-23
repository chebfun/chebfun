function fx = toValues(disc,f)
%TOVALUES  Convert chebfun to a COLLOC2 discretization.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% The easy part.
x = points(disc);
fx = f(x);

% Evaluate left- and right-sided limits at breaks.
n = disc.dimension;
csn = [0, cumsum(n)];
dxloc = csn(2:end-1);
fx(dxloc) = feval(f, x(dxloc), 'left');
fx(dxloc+1) = feval(f, x(dxloc), 'right');

end
