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

% Throw a warning if the result is not happy and we find roots in the domain:
if ( f.ishappy && ~g.ishappy )
    r = roots(f);
    if ( ~isempty(r) )
        warning(['Attempting to SQRT a CHEBTECH with one or more roots. ', ...
            'Result may be inaccurate.']);
    end
end

end