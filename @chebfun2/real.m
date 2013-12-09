function f = real(f)
%REAL  real part of a chebfun2.

if ( isempty( f ) )   % check for empty chebfun2.
    f = chebfun2;
    return
end

op = @(x,y) real( feval( f, x, y ) );  % Resample.
f = chebfun2( op, f.domain );           % Call constructor.

end