function F = diag(f, d)
%DIAG      Pointwise multiplication operator.
%   This function is deprecated. Use OPERATORBLOCK.MULT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = linop( operatorBlock.mult(f, double(d)) );

end
