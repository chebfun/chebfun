function g = power(f, b)
% .^   CHEBTECH power.
%   F.^G returns a CHEBTECH F to the scalar power G, a scalar F to the CHEBTECH
%   power G, or a CHEBTECH F to the CHEBTECH power G. F and or G may be complex.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If G is a constant, cast it to a CHEBTECH:
if ( isnumeric(b) )
    b = f.make(@(x) 0*x+b, b, 1);
end

% If F is complex, we need to do extrapolation at the end points in order to
% avoid the wrong sign due to the rounding errors, as F is supposed to have
% vanishing values at the end points.

pref = f.techPref();
if ( ~isreal(f) )
    pref.extrapolate = 1;
end

% Simply call the compose function:
g = compose(f, @power, b, pref);

% Throw a warning if the result is not happy and we find roots in the domain:
if ( f.ishappy && ~g.ishappy )
    r = roots(f);
    if ( ~isempty(r) )
        warning(['Attempting to compute the POWER of a CHEBTECH with one ', ...
            'or more roots. Result may be inaccurate.']);
    end
end

end