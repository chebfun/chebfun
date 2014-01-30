function f = imag(f)
%IMAG   Imaginary part of a SINGFUN.
%   IMAG(F) is the imaginary part of F.
%
% See also REAL, ISREAL, CONJ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Compute the imaginary part of the values:
f.funPart = imag(f.funPart);

% TODO: Is this right?
f.deltaMag = imag(f.deltaMag);
f = simplifyDeltas(f);

end
