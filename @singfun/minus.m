function s = minus(f,g)
%-   Subtraction of two SINGFUN objects.
%   F - G subtracts G from F, where F and G are SINGFUN objects or scalars.
%
% See also PLUS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call PLUS():
s = plus(f, -g);

end
