function f = minus(f, g)
%-   Subtraction of two CHEBTECH objects.
%   F - G subtracts G from F, where F and G are CHEBTECH objects or scalars.
%
% See also PLUS, UMINUS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% f - g = f + (-g)
f = plus(f, uminus(g)); 

end
