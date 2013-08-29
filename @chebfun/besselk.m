function g = besselk(nu, k, f, scale)
%BESSELK    Modified Bessel function of second kind of a CHEBFUN.
%   K = BESSELK(NU, F) computes the modified Bessel function of second kind
%   H1_nu(F) or H2_nu(F) of the nonzero CHEBFUN F. If F passes through the
%   origin in its domain, then an error is returned. The order NU need not be an
%   integer, but must be real. The argument F can be complex
%
%   K = BESSELK(NU, F, SCALE) returns a scaled K_nu(F) specfied by SCALE:
%         0 - (default) is the same as BESSELK(NU, F),
%         1 - scales H2_nu(F) by exp(F)).
%
% See also AIRY, BESSELH, BESSELI, BESSELJ, BESSELY.
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
    error('CHEBFUN:besselk:zero', 'F has roots in its domain.');
end

% Compose:
g = compose(f, @(x) besselk(nu, x));

% Scale (as described in help documentation):
if ( scale == 1 )
    g = exp(f).*g;
end

end
