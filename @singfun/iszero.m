function out = iszero(f)
%ISZERO   True for zero SINGFUN objects.
%   ISZERO(F) returns logical TRUE if the smooth part of  F is zero and FALSE
%   otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

out = iszero(f.smoothPart);

end
