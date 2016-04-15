function f = cosh( f )
%COSH   Hyperbolic cosine of a CHEBFUN3T.
%   COSH(F) returns the hyperbolic cosine of F.

% Check for empty:
if ( isempty( f ) ) 
    return 
end 

op = @(x,y,z) cosh( feval(f, x, y, z ) ) ;  % Resample. 
f = chebfun3t( op, f.domain );          % Call constructor. 

end