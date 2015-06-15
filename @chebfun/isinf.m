function out = isinf(f)
%ISINF   Test if a CHEBFUN is infinite.
%   ISINF(F) returns TRUE if F has any infinite values and FALSE otherwise.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = ~isfinite(f);

end
