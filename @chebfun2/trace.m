function t = trace( f )
% TRACE integral of a CHEBFUN2 along its diagonal 
%
% TRACE(f) is the integral of function f(x,x).
% 
% See also DIAG.

% Copyright 2013 by The University of Oxford and The Chebfun2 Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun2 information. 

t = sum( diag( f ) ); 

end