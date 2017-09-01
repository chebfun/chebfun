function out = isempty(f)
%ISEMPTY   True for an empty CLASSICFUN.
%   ISEMPTY(F) returns TRUE if F is an empty CLASSICFUN and FALSE otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if the ONEFUN is empty:
out = (numel(f) <= 1) && isempty(f.onefun);

end
