function f = floor(f)
%FLOOR   Pointwise floor function of a CHEBTECH.
%   G = FLOOR(F) returns the CHEBTECH G such that G(X) = FLOOR(F(x)) for each x
%   in [-1, 1].
%
%   If F is complex, then the G = FLOOR(REAL(F)) + 1i*FLOOR(IMAG(F)).
%
%   Note that FLOOR() assumes the output G(X) is a constant. If it is not, then
%   garbage is returned with no warning.
%
% See also CEIL, ROUND, FIX.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Evaluate at the two end points, and an arbitrary interior point:
arbitraryPoint = 0.1273881594;
fx = feval(f, [-1 ; arbitraryPoint ; 1]);
% Take the mean:
meanfx = mean(fx, 1);
% Compute the floor:
f.coeffs = floor(meanfx);

end
