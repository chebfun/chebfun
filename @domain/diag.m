function F = diag(f,d)
%DIAG      Pointwise multiplication operator.
%   This function is deprecated. Use OPERATORBLOCK.MULT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

F = linop( operatorBlock.mult(f) );

end
