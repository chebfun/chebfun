function f = real(f)
%REAL   Real part of a DELTAFUN.
%   REAL(F) is the real part of F.
%
% See also IMAG, ISREAL, CONJ.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Compute the real part of the classical part of F.
f.funPart = real(f.funPart);

end
