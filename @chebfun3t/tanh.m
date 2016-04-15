function f = tanh(f)
%TANH   Hyperbolic tangent of a CHEBFUN3T.
%   TANH(F) returns the hyperbolic tangent of a CHEBFUN3T.

% Empty check: 
if ( isempty( f ) ) 
    return
end

op = @(x,y,z) tanh( feval( f, x, y, z ) );  % Resample.    
f = chebfun3t( op, f.domain );          % Call constructor.

end
