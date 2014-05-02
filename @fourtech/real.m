function f = real(f)
%REAL   Real part of a FOURTECH.
%   REAL(F) is the real part of F.
%
%   See also ISREAL, IMAG, CONJ.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Compute the real part of the values:
f.values = real(f.values);
f.vscale = max(abs(f.values), [], 1);

if ( ~any(f.values(:)) )
    % Input was imaginary, so output a zero FOURTECH:
    f = f.make(zeros(1, size(f.values, 2)), f.vscale, f.hscale);
    f.ishappy = 1;
else
    % Compute real part of the coefficients:
    f.coeffs = f.vals2coeffs(f.values);
end

end
