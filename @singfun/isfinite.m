function out = isfinite(f)
%ISFINITE   Test if a SINGFUN is bounded.
%   ISFINITE(F) returns FALSE if F has any non trivial EXPONENT values and TRUE otherwise.
%
% See also ISINF, ISNAN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if F has exponents
% [TODO]: should we use the following:
% out = any(f.exponents); OR
out = all(abs(f.exponents) < singfun.pref.singfun.eps);
end
