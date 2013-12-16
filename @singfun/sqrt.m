function f = sqrt(f)
%SQRT   Square root of a SINGFUN.
%   SQRT(F) returns the square root of a SINGFUN F. Note, it is assumed that the
%   only roots of F are located at the endpoints of F.domain.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Extract roots from the boundaries and incrememnt the exponents accordingly:
f = extractBoundaryRoots(f);

% Exponents are halved by sqrt:
f.exponents = f.exponents/2;

% Call SQRT of f.smoothPart (the output of which is expected to be smooth):
f.smoothPart = sqrt(f.smoothPart);

end