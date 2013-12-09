function f = tanh(f)
%TANH Hyperbolic tangent of a chebfun2.
%
% TANH(F) returns the hyperbolic tangent of a chebfun2.

if ( isempty( f ) )  % check for empty chebfun2.
    return;
end

op = @(x,y) tanh( feval( f, x, y ) );  % Resample.    
f = chebfun2( op, f.domain );          % Call constructor.

end