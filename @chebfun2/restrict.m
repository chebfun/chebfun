function f = restrict(f, dom)
% RESTRICT  Restrict the domain of a chebfun2.
%
% F = RESTRICT(F, DOM) approximates the chebfun2 on the domain DOM. 

% Copyright 2013 by The University of Oxford and The Chebfun2 Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun2 information. 

f.cols = restrict(f.cols, dom(3:4)); 
f.rows = restrict(f.rows, dom(1:2));
f.domain = dom;

f = simplify( f ); 

end