function f = simplify( f ) 
% Simplify a SEPARABLEAPPROX
% 
% F = SIMPLIFY( F ) compressed the representation of F to one that is
% numerically the same, but requires fewer parameters to store. Currently this
% simplifies the polynomial degree of F, but not the rank. 

% Copyright 2015 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

% Simplify the column and row slices. 
f.cols = simplify( f.cols, [], 'globaltol' ); 
f.rows = simplify( f.rows, [], 'globaltol' ); 

% Note that we do not simplify the rank here because that would require 
% calling the constructor and is expensive.  In almost all situations we 
% are happy with simplifying the columns and rows only. 

end 
