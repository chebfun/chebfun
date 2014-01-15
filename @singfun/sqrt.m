function f = sqrt(f)
%SQRT   Square root of a SINGFUN.
%   SQRT(F) returns the square root of a SINGFUN F.
%
%   This function assumes that the curve traced out by F in the complex plane
%   both (1) does not come too close to zero except at the domain boundaries 
%   +/- 1 and (2) does not cross over the branch cut in SQRT along the negative
%   real axis.  That is, F should not vanish at any point of (-1, 1), and the
%   imaginary part of F should not vanish at any point of (-1, 1) where the real
%   part of F is negative.  If any of these assumptions are violated, garbage
%   may be returned with no warning.
%
% See also POWER.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Extract roots from the boundaries and incrememnt the exponents accordingly:
f = extractBoundaryRoots(f);

% Exponents are halved by sqrt:
f.exponents = f.exponents/2;

% Call SQRT of f.smoothPart (the output of which is expected to be smooth):
f.smoothPart = sqrt(f.smoothPart);

end