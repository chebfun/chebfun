function t = trace( f )
% TRACE integral of a chebfun2 along its diagonal 
%
% TRACE(f) is the integral of function f(x,x).
% 
% See also DIAG.

t = sum( diag( f ) ); 

end