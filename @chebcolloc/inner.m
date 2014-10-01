function F = inner(disc, f)
%INNER   Inner product functional for CHEBCOLLOC.
%   INNER(DISC, F) returns a row vector. The dot product of this vector with a 
%   CHEBCOLLOC2 discretization vector V results in the inner product of F with 
%   the function associated with V (at fixed discretization size). 
%
%   NOTE: The domain of F should match that of DISC. They are NOT checked, for
%   performance reasons.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get nodes and Clenshaw-Curtis quadrature weights:
[x, w] = functionPoints(disc);

% Curtis-Clenshaw weights times function values:
F = w.*f(x.');    

end
