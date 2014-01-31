function out = issmooth(f)
%ISSMOOTH   True if a SINGFUN object F is smooth.
%   ISSMOOTH(F) returns TRUE if the SINGFUN objects F has zero EXPONENTS or if 
%   the SMOOTHPART is zero. The test is FALSE otherwise.
%
%   See also ISSING.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

out = all(f.exponents == 0) || iszero(f.smoothPart);

end