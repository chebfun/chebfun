function f = sin( f ) 
% SIN   Sine of a chebfun2.

op = @(x,y) sin( feval(f, x, y) );  % Resample. 
f = chebfun2( op, f.domain );          % Call constructor.

end