function f = cancelExponents(f)
%CANCELEXPONENTS  Cancel the exponents of a SINGFUN.
%   F = CANCELEXPONENTS(F) returns a SINGFUN whose smoothPart has no vanishing
%   value at the boundaries if the corresponding exponents are negative.
%   Boundary roots are extracted off from the smoothPart at the corresponding
%   boundary to cancel the negative exponents as far as possible.
%
% See also EXTRACTBOUNDARYROOTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab the exponents:
exps = get(f, 'exponents');

% Grab boundary values:
boundaryVals = [lval(f.smoothPart),  rval(f.smoothPart)];

% Tolerance:
tol = 100*eps*vscale(f.smoothPart);

% Boundaries with negative exponent and vanishing value:
ind = ( exps < 0 ) & ( abs(boundaryVals) < repmat(tol, 1, 2) );

% Extract boundary roots to cancel the exponents:
if ( any(ind) )
    numRoots = (-exps.*ind).';
    f = extractBoundaryRoots(f, numRoots);
end

end
