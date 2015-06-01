function h = ge(f, g)
%>=   Greater than or equal operator for CHEBFUN objects.
%   H = F >= G, where F and/or G are CHEBFUN objects, constructs a logical
%   CHEBFUN H which is true (i.e., takes the value 1) where F >= G, and false
%   (0) elsewhere.
%
% See also LE, GT, LT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

h = le(g, f);

end
