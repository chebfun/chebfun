function f = not(f)
%~   CHEBTECH logical NOT.
%   NOT(F) returns a CHEBTECH which evaluates to one at all points where F is
%   zero and one otherwise.  F cannot have any roots in its domain.  If F does
%   have roots, then NOT(F) will return garbage with no warning. F may be
%   complex.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% TODO:  Should we use a tolerance here instead of any()?
f.values = ~any(f.values, 1);
f.coeffs = f.values;
f.vscale = abs(f.values);

end
