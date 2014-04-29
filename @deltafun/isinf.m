function out = isinf(f)
%ISINF   Test if a DELTAFUN is infinite.
%   ISINF(F) returns TRUE if F has any delta function or if the smooth part
%   is infinite. It is FALSE otherwise.
%
% See also ISFINITE, ISNAN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Simply call isinf
out = ~isfinite(f);

end
