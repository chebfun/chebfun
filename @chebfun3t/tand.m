function f = tand( f )
%TAND  Tangent of a CHEBFUN3T (in degrees)

% Empty check: 
if ( isempty( f ) )    
    return
end

op = @(x,y,z) tand( feval( f, x, y, z ) );   % Resample.
f = chebfun3t( op, f.domain );           % Call constructor.

end
