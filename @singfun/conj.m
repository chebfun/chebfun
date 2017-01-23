function f = conj(f)
%CONJ   Complex conjugate of a SINGFUN.
%   CONJ(F) is the complex conjugate of F. For a complex F,
%   CONJ(F) = REAL(F) - 1i*IMAG(F).
%
% See also REAL, IMAG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty arguments:
if ( isempty(f) )
    return;
end

% Conjugate the smooth part of F:
f.smoothPart = conj(f.smoothPart);

% Return a SMOOTHFUN object if F is smooth:
if ( issmooth(f) )
    f = f.smoothPart;
end

end
