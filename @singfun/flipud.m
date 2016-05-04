function f = flipud(f)
%FLIPUD   Flip/reverse a SINGFUN object.
%   G = FLIPUD(F) returns G such that G(x) = F(-x) for all x in [-1,1].

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Flip the smooth part:
f.smoothPart = flipud(f.smoothPart);

% Since exponents of a SINGFUN are contained in a 1x2 row vector, FLIPUD is 
% translated into a FLIPLR:
f.exponents = fliplr(f.exponents);

end
