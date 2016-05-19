function F = ctranspose(F)
%'   Conjugate transpose of a CHEBFUN3V object

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Transpose and then conjugate: 
F = transpose(F); 
F = conj(F); 

end