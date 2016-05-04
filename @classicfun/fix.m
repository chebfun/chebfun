function f = fix(f)
%FIX   Round a CLASSICFUN pointwise toward zero.
%   G = FIX(F) returns the CLASSICFUN G such that G(x) = FIX(F(x)) for each x in
%   F.domain.
%
%   If F is complex, then the G = FIX(REAL(F)) + 1i*FIX(IMAG(F)).
%
%   Note that FIX() assumes the output G(X) is a constant. If it is not, then
%   garbage is returned with no warning.
%
% See also ROUND, CEIL, FLOOR.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f.onefun = fix(f.onefun);

end
