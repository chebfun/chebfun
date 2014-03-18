function f = conj(f)
%CONJ   Complex conjugate of a FOURIERTECH.
%   CONJ(F) is the complex conjugate of F. For a complex F,
%   CONJ(F) = REAL(F) - 1i*IMAG(F).
%
% See also REAL, IMAG.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% No need to conjugate a real function
if f.isReal
    return;
end

% Conjugate the values:
f.values = conj(f.values);

% This won't quite work for the coefficients:
% % Conjugate the coefficients:
% f.coeffs = conj(f.coeffs);
% So instead just recompute the coefficients for the conjugated values.
f.coeffs = f.vals2coeffs(f.values);

end
