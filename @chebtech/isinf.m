function out = isinf(f)
%ISINF   Test if a CHEBTECH is unbounded.
%   ISINF(F) returns TRUE if F has any infinite values and FALSE otherwise.
%
% See also ISFINITE, ISNAN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if any coefficients are infinite:
out = any(isinf(f.coeffs(:)));

end
