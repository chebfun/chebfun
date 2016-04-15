function f = sin( f ) 
%SIN   Sine of a CHEBFUN3T.
%   SIN(F) returns the sine of F.
%

% Check for empty:
if ( isempty( f ) )
    return
end 

op = @(x,y, z) sin( feval( f, x, y, z ) );  % Resample. 
f = chebfun3t(op, f.domain);           % Call constructor. 

end
