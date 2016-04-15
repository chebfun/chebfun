function f = cos( f ) 
%COS   Cosine of a CHEBFUN3T.
%   COS(F) returns the cosine of F.
%

% Check for empty:
if ( isempty( f ) )
    return
end 

op = @(x,y,z) cos( feval( f, x, y, z ) );  % Resample. 
f = chebfun3t(op, f.domain);           % Call constructor. 

end
