function f = fliplr(f)
%FLIPLR   Flip columns of an array-valued CHEBTECH object.
%   FLIPLR(F) flips the columns of an array-valued CHEBTECH in the left/right
%   direction. If F has only one column, then this function has no effect.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Flip the orders of the columns of the matrices storing the coefficients.
f.coeffs = fliplr(f.coeffs);
f.vscale = fliplr(f.vscale);
f.epslevel = fliplr(f.epslevel);

end
