function f = flipud(f)
%FLIPUD   Flip/reverse a SINGFUN object.
%   G = FLIPUD(F) returns G such that G(x) = F(-x) for all x in [-1,1].

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Flip the smooth part:
f.smoothPart = flipud(f.smoothPart);

% Since the following fields of a SINGFUN are 1x2 row vectors, FLIPUD is 
% translated into a FLIPLR:
f.exponents = fliplr(f.exponents);
f.singType = fliplr(f.singType);

end
