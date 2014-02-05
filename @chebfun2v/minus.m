function f = minus( f, g )
% - MINUS. Minus of two chebfun2v.  
%   F - G substracts the chebfun2v F from G componentwise. 
%
%   minus(F, G) is called for the syntax f - g.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

f = plus( f, uminus( g ) );

end