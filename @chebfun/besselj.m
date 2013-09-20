function g = besselj(nu, f, scale, pref)
%BESSELJ    Bessel function of first kind of a CHEBFUN.
%   J = BESSELJ(NU, F) returns J_nu(F), i.e., is the Bessel function of the
%   first kind, J_nu(Z) composed with the chebfun object F. The order NU need
%   not be an integer, but must be a real scalar. The CHEBFUN F can be complex.
%
%   J = BESSELJ(NU, F, SCALE) returns a scaled J_nu(F) specified by SCALE:
%         0 - (default) is the same as besselj(NU, F)
%         1 -  scales J_nu(F) by exp(-abs(imag(F)))
%
% See also AIRY, BESSELH, BESSLI, BESSELK, BESSELY.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin < 4 )
    pref = chebfun.pref();
end

if ( any(size(nu) > 1) || ~isreal(nu) )
    error('CHEBFUN:besselj:nu', ...
        'The first argument of besselj must be a real-valued scalar.');
end

% Singular part:
fnu = f.^nu;

% Smooth part (see below):
g = compose(f, @(x) h(nu, x), pref);

% Combine:
g = fnu.*g;

% Scale (as described in help documentation):
if ( nargin == 3 && scale == 1 )
    scl = exp(-abs(imag(f)));
    g = scl.*g;
end

end

function y = h(nu, x)
    % h(nu,x) = J(nu,x)/x^nu is smooth!
    y = besselj(nu, x) ./ (x.^nu);
    % Deal with origin:
    x0 = x == 0;
    if ( any(x0) )
        y(x0) = 2^(-nu) / gamma(nu + 1);
    end
end
