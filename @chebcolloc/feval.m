function E = feval(disc, location, direction)
%FEVAL   Evaluation functional for CHEBCOLLOC.
%   FEVAL(DISC, LOC, DIRN) returns a functional that evaluates the Chebyshev
%   polynomial represented by a COLLOC discretization at the given point LOC as
%   approached from the direction DIRN (either +1 or -1).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = disc.dimension;

% Find the collocation points and create an empty functional.
[x, ignored, v] = functionPoints(disc);
offset = cumsum([0 ; n(:)]);
N = offset(end);
E = zeros(1, N);

% Only one subinterval creates nonzero entries in E.
intnum = disc.whichInterval(location, direction);
active = offset(intnum) + (1:n(intnum));
E(1, active) = barymat(location, x(active), v(active));

end
