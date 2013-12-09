function f = tan( f )
% TAN   Tangent of a chebfun2.

if ( isempty( f ) ) % check for empty chebfun2.
    return
end

op = @(x,y) tan( feval( f, x, y ) ); % resample
f = chebfun2( op, f.domain );        % Call constructor.

end