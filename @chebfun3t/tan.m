function f = tan( f )
%TAN   Tangent of a CHEBFUN3T.

% Empty check: 
if ( isempty( f ) ) 
    return
end

op = @(x,y,z) tan( feval( f, x, y, z ) ); % Resample
f = chebfun3t( op, f.domain );        % Call constructor.

end
