function f = real(f)
%REAL   Real part of a TRIGTECH.
%   REAL(F) is the real part of F.
%
% See also ISREAL, IMAG, CONJ.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If f is real then there is nothing to do:
if ( isreal(f) )
    return;
end

% Compute the real part of the values and update vscale:
f.values = real(f.values);

if ( ~any(f.values(:)) )
    % Input was imaginary, so output a zero TRIGTECH:
    f = f.make(zeros(1, size(f.values, 2)));
    f.ishappy = 1;
else
    % Compute the coefficients.
    f.coeffs = f.vals2coeffs(f.values);
end

f.isReal = true(1, size(f, 2));

% Simplify the result in case the imaginary part was contributing more to
% the length.
f = simplify(f);

end
