function B = mpower(A,pow)
%^         Repeated application of a linop.
%   B = A^m, for positive integer m, returns a linop that is the mth power
%   of A. Note that boundary conditions are not carried into the result. 

%  Copyright 2017 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

B = linop( mpower@chebmatrix(A,pow) );

end
