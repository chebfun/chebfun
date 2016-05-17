function out = isinf(f)
%ISINF   Test if a SINGFUN is unbounded.
%   ISINF(F) returns TRUE if F has any negative EXPONENT values and FALSE
%   otherwise.
%
% See also ISFINITE, ISNAN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Simply call and negate isfinite:
out = ~isfinite(f);

end
