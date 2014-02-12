function f = imag(f)
%IMAG   Imaginary part of a CHEBTECH.
%   IMAG(F) is the imaginary part of F.
%
% See also REAL, ISREAL, CONJ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Compute the imaginary part of the values:
f.values = imag(f.values);
f.vscale = max(abs(f.values), [], 1);

if ( ~any(f.values(:)) )
    % Input was real, so output a zero CHEBTECH:
    f = f.make(zeros(1, size(f.values, 2)), 0, f.hscale);
    f.ishappy = 1;
else
    % Compute imaginary part of the coefficients:
    f.coeffs = imag(f.coeffs);
end

end
