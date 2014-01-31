function f = extractBoundaryRoots(f)
%EXTRACTBOUNDARYROOTS   Extract boundary roots of the smooth part of a SINGFUN 
%   and absorb them into its exponents.
%
%   F = EXTRACTBOUNDARYROOTS(F) returns a SINGFUN G whose SMOOTHPART is free of 
%   roots at the boundary points -1 and 1 and these roots are represented by 
%   proper change in its EXPONENTS.
%
% See also ROOTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Extract the roots from the boundaries:
[f.smoothPart, rootsL, rootsR] = extractBoundaryRoots(f.smoothPart);

% Increment the exponents:
f.exponents = f.exponents + [rootsL, rootsR];

end