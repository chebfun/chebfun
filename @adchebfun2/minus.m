function f = minus(f,g)
% -	  Minus 
% 
% F-G subtracts ADchebfun2 G from F, or a scalar from a chebfun2 if either
% F or G is a scalar.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

% f-g =  f + (-g) 
f = plus(f,uminus(g));

end