function f = real(f)
%REAL	Real part of a FUNCHEB1.
%   REAL(F) is the real part of F.
%
%   See also ISREAL, IMAG, CONJ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Compute the real part of the values:
f.values = real(f.values);

if ( ~any(f.values(:)) )
    
    % Input was imaginary, so output a zero FUNCHEB1.
    f = funcheb1(zeros(1, size(f.values, 2)), f.vscale, f.epslevel);
    
else
    
    % Compute real part of the coeffs
    f.coeffs = real(f.coeffs);

end

end
