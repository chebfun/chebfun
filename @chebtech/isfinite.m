function out = isfinite(f)
%ISFINITE   Test if a CHEBTECH is bounded.
%   ISFINITE(F) returns FALSE if F has any infinite values and TRUE otherwise.
%
% See also ISINF, ISNAN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if coefficients are finite:
out = all(isfinite(f.coeffs(:)));

end
