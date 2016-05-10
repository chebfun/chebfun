function g = ceil(f)
%CEIL   Pointwise ceiling of a CHEBFUN.
%   G = CEIL(F) returns the CHEBFUN G such that G(x) = CEIL(F(x)) for each x in
%   F.domain.
%
%   If F is complex, then the G = CEIL(REAL(F)) + 1i*CEIL(IMAG(F)).
%
% See also FLOOR, ROUND, FIX.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with unbounded functions:
if ( ~isfinite(f) )
    error('CHEBFUN:CHEBFUN:ceil:inf', ...
        'Ceil is not defined for functions which diverge to infinity.');
end

g = -floor(-f);

end
