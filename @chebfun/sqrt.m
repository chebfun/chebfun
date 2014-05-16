function f = sqrt(f, pref)
%SQRT   Square root of a CHEBFUN.
%   SQRT(F) returns the square root of a CHEBFUN F.
%
% See also POWER.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    pref = chebfunpref();
end

% Trivial case: (f is empty)
if ( isempty(f) )
    return
end

% Simply call POWER()
f = power(f, 0.5, pref);

end
