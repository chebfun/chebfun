function f = tanh(f)
%TANH   Hyperbolic tangent of a CHEBFUN3.
%   tanh(F) returns the hyperbolic tangent of a CHEBFUN3.

% Empty check: 
if ( isempty(f) ) 
    return
end

op = @(x,y,z) tanh(feval(f, x, y, z));  % Resample.    
f = chebfun3(op, f.domain, 'fiberDim', 3);          % Call constructor.

end