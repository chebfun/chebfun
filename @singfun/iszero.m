function out = iszero(f)
%ISZERO   True for zero SINGFUN objects.
%   ISZERO(F) returns logical TRUE if the smooth part of F is zero and FALSE
%   otherwise.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if the smooth part is zero:
out = iszero(f.smoothPart);

end
