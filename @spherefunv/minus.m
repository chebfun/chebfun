function f = minus( f, g )
% - MINUS. Minus of two SPHEREFUNV.  
%   F - G subtracts the SPHEREFUNV F from G componentwise.
%
%   minus(F, G) is called for the syntax f - g.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


f = plus( f, uminus( g ) );

end
