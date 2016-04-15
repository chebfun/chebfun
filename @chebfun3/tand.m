function f = tand(f)
%TAND  Tangent of a CHEBFUN3 (in degrees)

% Empty check: 
if ( isempty(f) )
    return
end

op = @(x,y,z) tand(feval(f, x, y, z));               % Resample.
f = chebfun3(op, f.domain, 'fiberDim', 3);           % Call constructor.

end