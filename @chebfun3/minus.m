function f = minus(f, g)
%-   Subtraction of two CHEBFUN3 objects.
%   MINUS(F, G) subtracts G from F, where F and G are CHEBFUN3 objects.
%
% See also CHEBFUN3/PLUS and CHEBFUN3/UMINUS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% f - g = f + (-g)
f = plus(f, uminus(g));

end