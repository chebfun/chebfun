function f = minus(f, g)
%-	Subtraction of two CLASSICFUN objects.
%   F - G subtracts G from F, where F and G may be CLASSICFUN objects or scalars.
%
%   If F and G are both CLASSICFUN objects, they are assumed to have the same domain.
%   The method gives no warning if their domains don't agree, but the output of
%   the method will be meaningless.
%
% See also PLUS, UMINUS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% To avoid duplication, we note that (f - g) = f + (-g):
f = plus(f, uminus(g)); 

end
