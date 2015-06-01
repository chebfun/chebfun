function out = isempty(f)
%ISEMPTY   True for an empty DELTAFUN.
%   ISEMPTY(F) returns TRUE if F is an empty DELTAFUN and FALSE otherwise.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if the delta part and the Chebfun is empty:
out = isempty(f.deltaLoc) && isempty(f.funPart);

end
