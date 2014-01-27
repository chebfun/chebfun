function F = mult(disc,f)
%MULT      Multiplication operator in COLLOC2.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

F = diag( toValues(disc,f) );

end
