function g = besselh(nu, k, f, scale)
%BESSELH    Bessel function of third kind (Hankel function) of a CHEBFUN.
%   H = BESSELH(NU, K, F), for K = 1 or 2, computes the Hankel function H1_nu(F)
%   or H2_nu(F) of the nonzero CHEBFUN F. If F passes through the origin in its
%   domain, then an error is returned.
%
%   H = BESSELH(NU, F) uses K = 1.
%
%   H = BESSELH(NU, K, F, SCALE) returns a scaled Hankel function specfied by
%   SCALE:
%         0 - (default) is the same as BESSELH(NU, K, F)
%         1 - returns the following depending on K
%     H = BESSELH(NU, 1, F, 1) scales H1_nu(F) by exp(-i*F))).
%     H = BESSELH(NU, 2, F, 1) scales H2_nu(F) by exp(+i*F))).
%
%  The relationship between the Hankel and Bessel functions is:
%  
%         besselh(nu,1,z) = besselj(nu,z) + i*bessely(nu,z)
%         besselh(nu,2,z) = besselj(nu,z) - i*bessely(nu,z)
%
% See also AIRY, BESSELI, BESSELJ, BESSELK, BESSELY.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Parse inputs:
if ( nargin < 4 )
    scale = 1;
end
if ( isa(k, 'chebfun') )
    if ( nargin == 3 )
        scale = f;
    end
    f = k;
    k = 0;
end

% Check for roots:
r = roots(f, 'nojumps', 'nozerofun');
if ( numel(r) > 0 )
    error('CHEBFUN:besselh:zero', 'F has roots in its domain.');
end

% Compose:
g = compose(f, @(x) besselh(nu, k, x));

% Scale (as described in help documentation):
if ( scale == 1 )
    sgn = (-1)^k; % -1 for k = 1, 1 for k = 2.
    scl = exp(sgn*1i*f);
    g = scl.*g;
end

end
