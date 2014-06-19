function f = logical(f)
%LOGICAL   SINGFUN logical.
%   LOGICAL(F) returns a SMOOTHFUN which evaluates to one at all points where F
%   is nonzero and zero otherwise.  F cannot have any roots in its domain.  If
%   F does have roots, then LOGICAL(F) will return garbage with no warning. F
%   may be complex.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = logical(f.smoothPart);

end
