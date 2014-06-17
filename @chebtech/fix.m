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
% See http://www.chebfun.org/ for Chebfun information.

% Evaluate at the two end points, and an arbitrary interior point:
arbitraryPoint = 0.1273881594;
fx = feval(f, [-1 ; arbitraryPoint ; 1]);
% Take the mean:
meanfx = mean(fx, 1);
% Compute the fix:
f.coeffs = fix(meanfx);
f.vscale = abs(f.coeffs);
f.epslevel = 0*f.vscale + eps;

end
