function f = logical(f)
%LOGICAL   CLASSICFUN logical.
%   LOGICAL(F) returns a CLASSICFUN which evaluates to one at all points where F is
%   nonzero and zero otherwise.  F cannot have any roots in its domain.  If F
%   does have roots, then LOGICAL(F) will return garbage with no warning. F may
%   be complex.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f.onefun = logical(f.onefun);

end
