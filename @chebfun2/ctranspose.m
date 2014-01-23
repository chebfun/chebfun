function f = ctranspose( f )
%'	  Complex conjugate transpose of a chebfun2.
% 
% F' is the complex conjugate transpose of F.
%
% G = CTRANSPOSE(F) is called for the syntax F'.  

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Call transpose:
f = transpose( f ); 

% Take conj part:
f.cols = conj( f.cols ); 
f.rows = conj( f.rows ); 

end