function out = isempty(f)
%ISEMPTY   True for an empty CHEBTECH.
%   ISEMPTY(F) returns TRUE if F is an empty CHEBTECH and FALSE otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if the values are empty:
out = (numel(f) <= 1) && isempty(f.coeffs);

end
