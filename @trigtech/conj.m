function f = conj(f)
%CONJ   Complex conjugate of a TRIGTECH.
%   CONJ(F) is the complex conjugate of F. For a complex F,
%   CONJ(F) = REAL(F) - 1i*IMAG(F).
%
% See also REAL, IMAG.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% No need to conjugate a real function:
id = ~f.isReal;
if all(~id)
    return;
end

% Conjugate the values:
f.values(:, id) = conj(f.values(:, id));

% Could just recompute the coefficients for the conjugated values.
% f.coeffs = f.vals2coeffs(f.values);
% But this exploits the properties of the interpolant in complex exponential
% form:
f.coeffs(:, id) = flipud(conj(f.coeffs(:, id)));

end
