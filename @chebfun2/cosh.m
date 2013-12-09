function f = cosh( f )
%COSH Hyperbolic cosine of a chebfun2.
%
%  COSH(F) returns the hyperbolic cosine of F. 

if ( isempty( f ) ) % check for empty chebfun2.
    return; 
end 

op = @(x,y) cosh( feval(f, x, y ) ) ;  % Resample. 
f = chebfun2( op, f.domain );          % Call constructor. 

end