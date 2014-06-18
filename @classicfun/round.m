function f = round(f)
%ROUND   Round a CLASSICFUN pointwise to nearest integer.
%   G = ROUND(F) returns the CLASSICFUN G such that G(x) = ROUND(F(x)) for each x
%   in F.domain.
%
%   If F is complex, then the G = ROUND(REAL(F)) + 1i*ROUND(IMAG(F)).
%
%   Note that ROUND() assumes the output G(X) is a constant. If it is not, then
%   garbage is returned with no warning.
%
% See also CEIL, FLOOR, FIX.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% ROUND() the ONEFUN:
f.onefun = round(f.onefun);

end
