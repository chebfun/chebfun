function F = ctranspose( F )
% '   Conjugate transpose of a CHEBFUN2V

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

% Transpose and then conjugate: 
F = transpose( F ); 
F = conj( F ); 

end