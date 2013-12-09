function f = exp(f) 
% EXP  Exponential of a chebfun2

op = @(x,y) exp( feval(f, x, y) ); % resample.
f = chebfun2( op, f.domain );      % Call constructor.

end