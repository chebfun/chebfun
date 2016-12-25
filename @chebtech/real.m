function f = real(f)
%REAL   Real part of a CHEBTECH.
%   REAL(F) is the real part of F.
%
% See also ISREAL, IMAG, CONJ.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Compute the real part of the coefficients:
f.coeffs = real(f.coeffs);

if ( ~any(f.coeffs(:)) )
    % Input was imaginary, so output a zero CHEBTECH:
    f = f.make(zeros(1, size(f.coeffs, 2)));
    f.ishappy = 1;
end

end
