function f = fix(f)
%FIX   Round a CHEBTECH pointwise toward zero.
%   G = FIX(F) returns the CHEBTECH G such that G(x) = FIX(F(x)) for each x in
%   [-1, 1].
%
%   If F is complex, then the G = FIX(REAL(F)) + 1i*FIX(IMAG(F)).
%
%   Note that FIX() assumes the output G(X) is a constant. If it is not, then
%   garbage is returned with no warning.
%
% See also ROUND, CEIL, FLOOR.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

arbitraryPoint = 0.1273881594;
f.coeffs = fix(feval(f, arbitraryPoint));
f.vscale = abs(f.coeffs);

end
