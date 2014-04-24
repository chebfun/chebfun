function F = mult(disc, f)
%MULT   Multiplication operator in COLLOC.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = diag( toValues(disc, f) );

end
