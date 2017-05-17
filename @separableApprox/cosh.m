function f = cosh( f )
%COSH   Hyperbolic cosine of a SEPARABLEAPPROX.
%   COSH(F) returns the hyperbolic cosine of F.
% 
% See also SINH, COS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty( f ) ) 
    return 
end 

f = compose( f, @cosh ); 

end
