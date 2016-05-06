function f = fliplr(f)
%FLIPLR   Flip columns of an array-valued TRIGTECH object.
%   FLIPLR(F) flips the columns of an array-valued TRIGTECH in the left/right
%   direction. If F has only one column, then this function has no effect.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Flip the orders of the columns of the matrices storing the values and
% coefficients.
f.values = fliplr(f.values);
f.coeffs = fliplr(f.coeffs);

end
