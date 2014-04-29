function [f, rootsLeft, rootsRight] = extractBoundaryRoots(f, numRoots)
%EXTRACTBOUNDARYROOTS   Extract boundary roots of the smooth part of a SINGFUN 
%   and absorb them into its exponents.
%
%   [G, ROOTSLEFT, ROOTSRIGHT] = EXTRACTBOUNDARYROOTS(F) returns a SINGFUN G
%   which is free of roots at the boundary points -1 and 1. The multiplicity of
%   the boundary roots extracted at -1 and 1 are ROOTSLEFT and ROOTRIGHT, 
%   respectively.
%
%   [G, ROOTSLEFT, ROOTSRIGHT] = EXTRACTBOUNDARYROOTS(F, NUMROOTS) returns a 
%   SINGFUN G whose SMOOTHPART has been peeled off roots from the boundaries 
%   with the multiplicities specified by ROOTSLEFT and ROOTSRIGHT for the left 
%   boundary and the right boundary, respectively. NUMROOTS is a 2x1 vector 
%   specifying the multiplicities of the boundary roots that 
%   EXTRACTBOUNDARYROOTS aims to extract. The first and the second entry of 
%   NUMROOTS correspond to the left and the right boundary, respectively. Note 
%   that there is no warning thrown, when ROOTSLEFT and ROOTSRIGHT are not 
%   consistent with NUMROOTS.
%
% See also ROOTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    % Extract the roots from the boundaries
    [f.smoothPart, rootsLeft, rootsRight] = ...
        extractBoundaryRoots(f.smoothPart);
else
    % Extract the roots from the boundaries
    [f.smoothPart, rootsLeft, rootsRight] = ...
        extractBoundaryRoots(f.smoothPart, numRoots);
end

% Increment the exponents:
f.exponents = f.exponents + [rootsLeft, rootsRight];

end
