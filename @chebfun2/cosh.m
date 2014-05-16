function f = cosh( f )
%COSH   Hyperbolic cosine of a CHEBFUN2.
%   COSH(F) returns the hyperbolic cosine of F.
% 
% See also SINH, COS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty( f ) ) 
    return 
end 

op = @(x,y) cosh( feval(f, x, y ) ) ;  % Resample. 
f = chebfun2( op, f.domain );          % Call constructor. 

end