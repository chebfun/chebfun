function f = not(f)
%~   SINGFUN logical NOT.
%   NOT(F) returns a SMOOTHFUN which evaluates to one at all points where F is
%   zero and one otherwise.  F cannot have any roots in its domain.  If F does
%   have roots, then NOT(F) will return garbage with no warning. F may be
%   complex.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = ~f.smoothPart;

end
