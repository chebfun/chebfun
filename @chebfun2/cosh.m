function f = cosh( f )
%COSH   Hyperbolic cosine of a CHEBFUN2.
%   COSH(F) returns the hyperbolic cosine of F.
% 
% See also SINH, COS.

% Check for empty:
if ( isempty( f ) ) 
    return 
end 

op = @(x,y) cosh( feval(f, x, y ) ) ;  % Resample. 
f = chebfun2( op, f.domain );          % Call constructor. 

end