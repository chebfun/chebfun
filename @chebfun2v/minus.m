function f = minus( f, g )
% - MINUS. Minus of two chebfun2v.  
%
% f - g substracts the chebfun2v f from g componentwise. 
% minus(f,g) is called for the syntax f - g. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

f = plus( f, uminus( g ) );

end