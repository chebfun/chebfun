function f = imag(f)
%IMAG   Imaginary part of a TRIGTECH.
%   IMAG(F) is the imaginary part of F.
%
% See also REAL, ISREAL, CONJ.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isreal(f) )
    % Input was real, so output a zero TRIGTECH:
    f = f.make(zeros(1, size(f.values, 2)));
    f.ishappy = 1;
else
    % Compute the imaginary part of the values and update vscale:
    f.values = imag(f.values);
    % Compute imaginary part of the coefficients:
    f.coeffs = f.vals2coeffs(f.values);
    f.isReal = true(1,size(f.values,2));
end

end
