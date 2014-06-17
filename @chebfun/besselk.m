function F = besselk(nu, F, scale, pref)
%BESSELK   Modified Bessel function of second kind of a CHEBFUN.
%   K = BESSELK(NU, F) computes the modified Bessel function of second kind
%   K_NU(F) of the nonzero CHEBFUN F. If F passes through the origin in its
%   domain, then an error is returned. The order NU need not be an integer but
%   must be real. The argument F can be complex.
%
%   K = BESSELK(NU, F, SCALE) returns a scaled K_NU(F) specified by SCALE:
%         0 - (default) is the same as BESSELK(NU, F),
%         1 - scales K_NU(F) by exp(F)).
%
%   K = BESSELK(NU, F, SCALE, PREF) uses the CHEBFUNPREF object PREF when
%   building the CHEBFUN K.
%
% See also AIRY, BESSELH, BESSELI, BESSELJ, BESSELY.
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
    F(k) = columnBesselk(nu, F(k), scale, pref);
end

end

function g = columnBesselk(nu, f, scale, pref)

% Check for roots:
r = roots(f, 'nojump', 'nozerofun');
if ( numel(r) > 0 )
    error('CHEBFUN:CHEBFUN:besselk:zero', 'F has roots in its domain.');
end

% Compose:
g = compose(f, @(x) besselk(nu, x), pref);

% Scale (as described in help documentation):
if ( scale == 1 )
    scl = exp(f);
    g = scl.*g;
end

end
