function f = sin(f)
%SIN   Sine of a CHEBFUN3.
%   SIN(F) returns the sine of F.

% Check for empty:
if ( isempty(f) )
    return
end 

op = @(x,y, z) sin(feval(f, x, y, z));  % Resample. 
f = chebfun3(op, f.domain, 'fiberDim', 3);           % Call constructor. 

end