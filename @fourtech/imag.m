function f = imag(f)
%IMAG   Imaginary part of a FOURIERTECH.
%   IMAG(F) is the imaginary part of F.
%
%   See also REAL, ISREAL, CONJ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

if ( f.isReal )
    % Input was real, so output a zero FOURIERTECH:
    f = f.make(zeros(1, size(f.values, 2)), f.vscale, f.hscale);
    f.ishappy = 1;
else
    % Compute the imaginary part of the values:
    f.values = imag(f.values);
    % Compute imaginary part of the coefficients:
    f.coeffs = f.vals2coeffs(f.values);
    f.vscale = max(abs(f.values), [], 1);
end
f.isReal = true;

end
