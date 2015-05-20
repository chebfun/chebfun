function f = conj(f)
%CONJ   Complex conjugate of a DELTAFUN.
%   CONJ(F) is the complex conjugate of F. The FUNPART is conjugated as well as
%   the matrix of delta function magnitude.
%
% See also REAL, IMAG.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Conjugate the classical part:
f.funPart = conj(f.funPart);
f.deltaMag = conj(f.deltaMag);

end
