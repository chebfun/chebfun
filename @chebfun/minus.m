function f = minus(f, g)
%-	  Minus.
%   F - G subtracts chebfuns F and G, or a scalar from a chebfun if either F or
%   G is a scalar.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = plus(f, uminus(g));

end