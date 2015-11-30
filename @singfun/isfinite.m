function out = isfinite(f)
%ISFINITE   Test if a SINGFUN is bounded.
%   ISFINITE(F) returns FALSE if F has any non trivial EXPONENT values and
%   TRUE otherwise.
%
% See also ISINF, ISNAN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

tol = chebfunpref().blowupPrefs.exponentTol;

% Check if F has exponents and the smooth part is finite:
out = all(f.exponents > -tol) && isfinite(f.smoothPart);

end
