function F = mult(disc,f)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
F = diag( toValues(disc,f) );

end
