function f = extractBoundaryRoots(f)
%EXTRACTBOUNDARYROOTS   Extract roots of the SMOOTHPART of F at the boundary 
%   points -1 and 1 and absorb them to EXPONENTS.
%
%   G = EXTRACTBOUNDARYROOTS(F) returns a SINGFUN G, whose SMOOTHPART is free of
%   roots at the boundary points -1 and 1. The roots are absorbed into 
%   EXPONENTS.
%
% See also ROOTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[f.smoothPart, rootsLeft, rootsRight] = extractBoundaryRoots(f.smoothPart);
f.exponents = [f.exponents(1)-rootsLeft f.exponents(2)-rootsRight];

end