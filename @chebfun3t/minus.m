function f = minus(f, g)
%-   Subtraction of two CHEBFUN3T objects.
%   F - G subtracts G from F, where F and G are CHEBFUN3T objects or 
%   scalars.
%
% See also CHEBFUN3T/PLUS and CHEBFUN3T/UMINUS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% f - g = f + (-g)
f = plus(f, uminus(g)); 

end