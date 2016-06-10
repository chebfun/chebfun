function f = minus(f, g)
% - MINUS. Minus of two CHEBFUN3V.  
%   F - G subtracts the CHEBFUN3V F from G componentwise.
%
%   minus(F, G) is called for the syntax f - g.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = plus(f, uminus(g));

end