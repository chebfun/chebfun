function f = minus(f, g)
%-   CHEBFUN minus.
%   F - G subtracts CHEBFUNs F and G, or a scalar from a CHEBFUN if either F or
%   G is a scalar.
%
%   H = MINUS(F, G) is called for the syntax 'F - G'.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = plus(f, uminus(g));

end
