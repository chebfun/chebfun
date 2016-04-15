function f = sinh(f)
%SINH   Hyperbolic sine of a CHEBFUN3.
%   sinh(F) returns the hyperbolic cosine of F.

% Empty check: 
if ( isempty(f) )
    return
end 

op = @(x,y,z) sinh(feval(f, x, y, z));  % Resample.
f = chebfun3(op, f.domain, 'fiberDim', 3);           % Call constructor.

end