function f = tand(f)
% TAND   Tangent of a chebfun2 (in degrees)

if ( isempty(f) )    % check for empty chebfun2.
    return
end

op = @(x,y) tand( feval( f, x, y ) );  % Resample.
f = chebfun2( op, f.domain );           % Call constructor.

end