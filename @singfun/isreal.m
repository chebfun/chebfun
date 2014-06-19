function out = isreal(f)
%ISREAL   True for real SINGFUN.
%   ISREAL(F) returns TRUE if the smooth part of F is real and FALSE otherwise.
%
% See also REAL, IMAG.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if the smooth part is real:
out = isreal(f.smoothPart);

end
