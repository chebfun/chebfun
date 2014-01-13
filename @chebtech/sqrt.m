function g = sqrt(f)
%SQRT   Square root of a CHEBTECH.
%   SQRT(F) returns the square root of a CHEBTECH F. Note, it is assumed that F
%   is non-zero on its domain.
%
% See also POWER.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If F is complex, we need to do extrapolation at the end points in order to
% avoid the wrong sign due to the rounding errors, as F is supposed to have
% vanishing values at the end points.

pref = f.techPref();
if ( ~isreal(f) )
    pref.extrapolate = 1;
end

% Simply call the compose function:
g = compose(f, @sqrt, [], pref);

end