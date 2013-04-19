function f = fliplr(f)
%FLIPLR   Flip columns of a vectorised CHEBTECH object.
%   FLIPLR(F) flips the columns of a vectorised CHEBTECH in the left/right
%   direction. If F has only one column, then this function has no effect.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Flip the orders of the columns of the matrices storing the values and
% coefficients.
f.values = fliplr(f.values);
f.coeffs = fliplr(f.coeffs);

end
