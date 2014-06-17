function g = besseli(nu, f, varargin)
%BESSELI    Modified Bessel function of first kind of a CHEBFUN.
%   I = BESSELI(NU, F) returns I_NU(F), i.e., is the modified Bessel function of
%   the first kind, I_NU(Z) composed with the CHEBFUN object F. The order NU
%   need not be an integer, but must be a real scalar. The CHEBFUN F can be
%   complex.
%
%   I = besseli(NU, Z, SCALE) returns a scaled I_NU(Z) specified by SCALE:
%        0 - (default) is the same as besseli(NU, Z)
%        1 -  scales I_NU(Z) by exp(-abs(real(Z)))
%
%   The relationship between I_NU(Z) and J_NU(Z) = BESSELJ(NU, Z) is
%
%        I_NU(Z) = 1i^NU * J_NU(1i*Z).
%
% See also AIRY, BESSELH, BESSLJ, BESSELK, BESSELY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

scl = 1i.^nu;
g = scl*besselj(nu, 1i*f, varargin{:});

end
