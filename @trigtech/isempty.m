function out = isempty(f)
%ISEMPTY   True for an empty TRIGTECH.
%   ISEMPTY(F) returns TRUE if F is an empty TRIGTECH and FALSE otherwise.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if the values are empty:
out = (numel(f) <= 1) && isempty(f.values);

end
