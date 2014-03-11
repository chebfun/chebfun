function fx = toValues(disc,f)
%TOVALUES  Convert CHEBFUN to a COLLOC2 discretization.
%   TOVALUES(DISC,F) converts a chebfun F to values at 2nd kind points in
%   the DISC.DOMAIN.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% The easy part.
x = functionPoints(disc);
fx = f(x);

% Evaluate left- and right-sided limits at breaks, which are part of the
% discretization.
n = disc.dimension;
csn = [0, cumsum(n)];   % offsets into breakpoints
dxloc = csn(2:end-1);
fx(dxloc) = feval(f, x(dxloc), 'left');
fx(dxloc+1) = feval(f, x(dxloc), 'right');

end
