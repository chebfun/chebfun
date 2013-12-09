function f = conj(f)
%CONJ   Complex conjugate of a DELTAFUN.
%   CONJ(F) is the complex conjugate of F.
%
% See also REAL, IMAG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Conjugate the classical part:
f.funPart = conj(f.funPart);
% Invert the isConj field in delta functions:
f.delta.isConj = ~f.delta.isConj;
end
