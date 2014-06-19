function f = round(f)
%ROUND   Pointwise round function of a CHEBTECH.
%   G = FLOOR(F) returns the CHEBTECH G such that G(X) = FLOOR(F(x)) for each x
%   in [-1, 1].
%
%   If F is complex, then the G = ROUND(REAL(F)) + 1i*ROUND(IMAG(F)).
%
%   Note that ROUND() assumes the output G(X) is a constant. If it is not, then
%   garbage is returned with no warning.
%
% See also CEIL, FLOOR, FIX.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Evaluate at the two end points, and an arbitrary interior point:
arbitraryPoint = 0.1273881594;
fx = feval(f, [-1 ; arbitraryPoint ; 1]);
% Take the mean:
meanfx = mean(fx, 1);
% Compute the round:
f.coeffs = round(meanfx);
f.vscale = abs(f.coeffs);
f.epslevel = 0*f.vscale + eps;

end
