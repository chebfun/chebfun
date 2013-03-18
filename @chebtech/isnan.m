function out = isnan(f)
%ISNAN  Test if a CHEBTECH is has any NaN values.
%   ISNAN(F) returns true if F has any NaN values and false otherwise.
%
% See also ISFINITE, ISINF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if any values are NaN:
out = any(isnan(f.values(:)));

end
