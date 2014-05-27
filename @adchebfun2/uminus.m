function f = uminus(f)
% -	  Unary minus.
% 
% -F negates the ADchebfun2 F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

f.chebfun2 = uminus(f.chebfun2);

f.der = -f.der;
end