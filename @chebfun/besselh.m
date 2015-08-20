function F = besselh(nu, k, F, scale, pref)
%BESSELH   Bessel function of third kind (Hankel function) of a CHEBFUN.
%   H = BESSELH(NU, K, F), for K = 1 or 2, computes the Hankel function H1_NU(F)
%   or H2_NU(F) of the nonzero CHEBFUN F. If F passes through the origin in its
%   domain, then an error is returned.  The CHEBFUN F may be complex.
%
%   H = BESSELH(NU, F) uses K = 1.
%
%   H = BESSELH(NU, K, F, SCALE) returns a scaled Hankel function specified by
%   SCALE:
%         0 - (default) is the same as BESSELH(NU, K, F)
%         1 - returns the following depending on K
%     H = BESSELH(NU, 1, F, 1) scales H1_NU(F) by exp(-i*F))).
%     H = BESSELH(NU, 2, F, 1) scales H2_NU(F) by exp(+i*F))).
%
%   H = BESSELH(NU, K, F, SCALE, PREF) uses the CHEBFUNPREF object PREF when
%   building the CHEBFUN H.
%
%  The relationship between the Hankel and Bessel functions is:
%  
%         besselh(nu,1,z) = besselj(nu,z) + 1i*bessely(nu,z)
%         besselh(nu,2,z) = besselj(nu,z) - 1i*bessely(nu,z)
%
% See also AIRY, BESSELI, BESSELJ, BESSELK, BESSELY.
%
% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin < 5 )
    pref = chebfunpref();
end
if ( nargin < 4 )
    scale = 0;
end

if ( isa(k, 'chebfun') )
    if ( nargin == 3 )
        scale = F;
    end
    F = k;
    k = 1;
end

% Loop over the columns:
for j = 1:numel(F)
    F(j) = columnBesselh(nu, k, F(j), scale, pref);
end

end

function g = columnBesselh(nu, k, f, scale, pref)

% Check for roots:
r = roots(f, 'nojump', 'nozerofun');
if ( numel(r) > 0 )
    error('CHEBFUN:CHEBFUN:besselh:zero', 'F has roots in its domain.');
end

% Compose:
g = compose(f, @(x) besselh(nu, k, x), pref);

% Scale (as described in help documentation):
if ( scale == 1 )
    sgn = (-1)^k; % -1 for k = 1, 1 for k = 2.
    scl = exp(sgn*1i*f);
    g = scl.*g;
end

end
