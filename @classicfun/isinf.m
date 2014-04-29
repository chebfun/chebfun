function out = isinf(f)
%ISINF   Test if a CLASSICFUN is unbounded.
%   ISINF(F) returns TRUE if F takes an infinite value and FALSE otherwise.
%
% See also ISFINITE, ISNAN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if the ONEFUN of f is infinite:
out = isinf(f.onefun);

end
