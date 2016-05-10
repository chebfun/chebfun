function t = trace( f )
% TRACE integral of a SEPARABLEAPPROX along its diagonal 
%
% TRACE(f) is the integral of function f(x,x).
% 
% See also DIAG.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

t = sum( diag( f ) ); 

end
