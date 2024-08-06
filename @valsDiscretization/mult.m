function F = mult(disc, f)
%MULT   Multiplication operator in VALSDISCRETIZATION.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% F = diag( toValues(disc, f) );

v = toValues(disc, f);
N = length(v);
F = spdiags(v, 0, N, N);

end
