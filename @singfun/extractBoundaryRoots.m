function f = extractBoundaryRoots(f)
% Extract roots from ends of singfuns.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Extract the roots from the boundaries
[f.smoothPart, rootsL, rootsR] = extractBoundaryRoots(f.smoothPart);

% Increment the exponents:
f.exponents = f.exponents + [rootsL, rootsR];

end