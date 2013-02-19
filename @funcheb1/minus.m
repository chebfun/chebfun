function f = minus(f, g)
%-	Subtraction of two FUNCHEB1 objects.
%   F - G subtracts G from F, where F and G may be FUNCHEB1 objects or scalars.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% f1 - f2 = f1 + (-f2)
f = plus(f, uminus(g)); 

end

