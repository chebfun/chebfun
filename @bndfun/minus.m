function f = minus(f, g)
%-	Subtraction of two BNDFUN objects.
%   F - G subtracts G from F, where F and G may be BNDFUN objects or
%   scalars.
%
% See also PLUS, UMINUS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% f1 - f2 = f1 + (-f2)
f = plus(f, uminus(g)); 

end

