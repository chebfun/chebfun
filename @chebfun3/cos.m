function f = cos(f)
%COS   Cosine of a CHEBFUN3.
%   COS(F) returns the cosine of F.

% Check for empty:
if ( isempty(f) )
    return
end 

op = @(x,y, z) cos(feval(f, x, y, z));  % Resample. 
f = chebfun3(op, f.domain, 'fiberDim', 3);          % Call constructor. 

end