function f = sinh(f)
%SINH Hyperbolic sine of a chebfun2.

if ( isempty( f ) )     % check for empty chebfun2.
    return;
end 

op = @(x,y) sinh( feval( f, x, y ) );  % Resample.
f = chebfun2( op, f.domain );          % Call constructor. 

end