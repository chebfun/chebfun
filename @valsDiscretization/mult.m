function F = mult(disc, f)
%MULT   Multiplication operator in VALSDISCRETIZATION.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = diag( toValues(disc, f) );

end
