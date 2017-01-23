function s = minus(f,g)
%-   Subtraction of two DELTAFUN objects.
%   F - G subtracts G from F, where F and G are DELTAFUN objects or scalars.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

s = plus(f, -g);

end
