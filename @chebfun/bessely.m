function g = bessely(nu, f, scale, pref)
%BESSELY   Bessel function of second kind of a CHEBFUN.
%   Y = BESSELY(NU, F) computes the Bessel function of the second kind Y_NU(F)
%   of the nonzero CHEBFUN F. The order NU need not be an integer but must be
%   real. The argument F can be complex but must not pass through the origin.
%   The result is real where F is positive.
%
%   Y = BESSELY(NU, F, SCALE) returns a scaled Y_NU(F) specified by SCALE:
%        0 - (default) is the same as BESSELY(NU, F)
%        1 -  scales Y_NU(F) by exp(-abs(imag(F)))
%
%   Y = BESSELY(NU, F, SCALE, PREF) uses the preference structure PREF when
%   building the CHEBFUN Y.
%
% See also AIRY, BESSELH, BESSELI, BESSELJ, BESSELK.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin < 4 )
    pref = chebfun.pref();
end

% Check for roots:
r = roots(f, 'nojump', 'nozerofun');
if ( numel(r) > 0 )
    error('CHEBFUN:bessely:zero', 'F has roots in its domain.');
end

% Compose:
g = compose(f, @(x) bessely(nu, x), pref);

% Scale (as described in help documentation):
if ( (nargin >= 3) && (scale == 1) )
    scl = exp(-abs(imag(f)));
    g = scl.*g;
end

end
