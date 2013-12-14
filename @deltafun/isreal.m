function out = isreal(f)
%ISREAL   True for real SINGFUN.
%   ISREAL(F) returns TRUE if the smooth part of F is real and FALSE otherwise.
%
%   See also REAL, IMAG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Check if the funPart is real
out = isreal(f.funPart);
%[TODO]: What to do with deltafunctions

end
