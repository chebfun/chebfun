function F = inner(disc, f)
%INNER     Inner product functional for COLLOC2.

% INNER(DISC,F) returns a row vector. The dot product of this vector with a COLLOC2
% discretization vector V results in the inner product of F with the function
% associated with V (at fixed discretization size). 
%
% NOTE: The domain of F should match that of DISC. They are NOT checked, for
% performance reasons.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

[x, w] = points(disc);
F = w.*f(x.');    % Curtis-Clenshaw weights times function values

end
