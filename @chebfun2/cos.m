function f = cos( f ) 
% COS   Cosine of a chebfun2.

op = @(x,y) cos( feval( f, x, y ) );  % Resample. 
f = chebfun2(op, f.domain);          % Call constructor. 

end