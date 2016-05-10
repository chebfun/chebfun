function out = isnan(f)
%ISNAN   Test if the FUNPART of a DELTAFUN has any NaN values.
%   ISNAN(F) returns TRUE if the FUNPART of F has any NaN values and 
%   FALSE otherwise.
%
% See also ISFINITE, ISINF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if any values are NaN:
out = isnan(f.funPart);

end
