function g = bessely(nu, f, scale)
%BESSELY    Bessel function of second kind of a CHEBFUN.
%   H = BESSELY(NU, F), for K = 1 or 2, computes the Bessel function of the
%   second kind Y_NU(F) of the nonzero CHEBFUN F. The order NU need not be an
%   integer, but must be real. The argument F can be complex, but must not pass
%   through the origin. The result is real where F is positive.
%
%   Y = bessely(NU, F, SCALE) returns a scaled Y_NU(F) specified by SCALE:
%        0 - (default) is the same as bessely(NU, F)
%        1 -  scales Y_NU(F) by exp(-abs(imag(F)))
%
% See also AIRY, BESSELH, BESSELI, BESSELJ, BESSELK.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Check for roots:
% [TODO]: Uncomment once merged to development.
% r = roots(f, 'nojumps', 'nozerofun'); 
r = roots(f);
if ( numel(r) > 0 )
    error('CHEBFUN:bessely:zero', 'F has roots in its domain.');
end

% Compose:
g = compose(f, @(x) bessely(nu, x));

% Scale (as described in help documentation):
if ( nargin == 3 && scale == 1 )
    scl = exp(-abs(imag(f)));
    g = scl.*g;
end

end
