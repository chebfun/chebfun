function out = isnan(f)
%ISNAN   Test if a CHEBTECH has any NaN values.
%   ISNAN(F) returns TRUE if F has any NaN values and FALSE otherwise.
%
% See also ISFINITE, ISINF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if any coefficients are NaN:
out = any(isnan(f.coeffs(:)));

end
