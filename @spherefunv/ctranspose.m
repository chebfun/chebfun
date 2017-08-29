function F = ctranspose(F)
% '   Conjugate transpose of a SPHEREFUNV.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = transpose( F ); 
F = conj( F ); 

end
