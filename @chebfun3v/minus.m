function f = minus(f, g)
% - MINUS. Minus of two CHEBFUN3V objects.
%   F - G subtracts the CHEBFUN3V object G from F.
%
%   MINUS(F, G) is called for the syntax F - G.
%
% See also CHEBFUN3V/UMINUS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = plus(f, uminus(g));

end