function out = isnan(f)
%ISNAN   Test if a BNDFUN is has any NaN values.
%   ISNAN(F) returns TRUE if F has any NaN values and FALSE otherwise.
%
% See also ISFINITE, ISINF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if the onefun of f has any NaN values:
out = isnan(f.onefun);

end
