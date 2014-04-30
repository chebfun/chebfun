function f = real(f)
%REAL   Real part of a CHEBTECH.
%   REAL(F) is the real part of F.
%
%   See also ISREAL, IMAG, CONJ.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Compute the real part of the values:
% f.values = real(f.values);

f.coeffs = real(f.coeffs);
if ( ~any(f.coeffs(:)) )
    % Input was imaginary, so output a zero CHEBTECH:
    f = f.make(zeros(1, size(f.coeffs, 2)), f.vscale, f.hscale);
    f.ishappy = 1;
else
    % Compute real part of the coefficients:
    f.coeffs = real(f.coeffs);
end
f.vscale = getvscl(f);
end
