function out = isfinite(f)
%ISFINITE   Test if a SINGFUN is bounded.
%   ISFINITE(F) returns FALSE if F has any non trivial EXPONENT values and TRUE otherwise.
%
% See also ISINF, ISNAN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if F has exponents and the smooth part is finite.
out = all(f.exponents > -singfun.pref.singfun.exponentTol) && ...
    isfinite(f.smoothPart);

end
