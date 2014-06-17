function F = bessely(nu, F, scale, pref)
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
%   Y = BESSELY(NU, F, SCALE, PREF) uses the CHEBFUNPREF object PREF when
%   building the CHEBFUN Y.
%
% See also AIRY, BESSELH, BESSELI, BESSELJ, BESSELK.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 4 )
    pref = chebfunpref();
end
if ( nargin < 3 )
    scale = 0;
end

% Loop over the columns:
for k = 1:numel(F)
    F(k) = columnBessely(nu, F(k), scale, pref);
end

end

function g = columnBessely(nu, f, scale, pref)

% Check for roots:
r = roots(f, 'nojump', 'nozerofun');
if ( numel(r) > 0 )
    error('CHEBFUN:CHEBFUN:bessely:zero', 'F has roots in its domain.');
end

% Compose:
g = compose(f, @(x) bessely(nu, x), pref);

% Scale (as described in help documentation):
if ( scale == 1 )
    scl = exp(-abs(imag(f)));
    g = scl.*g;
end

end
