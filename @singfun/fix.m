function g = fix(f)
%FIX   Round a SINGFUN pointwise toward zero.
%   G = FIX(F) returns the SMOOTHFUN G such that G(x) = FIX(F(x)) for each x in
%   [-1, 1].
%
%   If F is complex, then the G = FIX(REAL(F)) + 1i*FIX(IMAG(F)).
%
%   Note that FIX() assumes the output G(X) is a constant. If it is not, then
%   garbage is returned with no warning. Also note that FIX(F) makes sense
%   only when F is a finite SINGFUN, otherwise G(X) is garbage. 
%
% See also ROUND, CEIL, FLOOR.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

arbitraryPoint = 0.1273881594;
value = fix(feval(f, arbitraryPoint));
g = f.constructSmoothPart(@(x) value + 0*x, value, 1, []);

end