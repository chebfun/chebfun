function h = gt(f, g)
% GT (>) for CHEBFUN object.
%   H = F > G, where F and/or G are CHEBFUN objects, constructs a logical
%   CHEBFUN H which is true (i.e., takes the value 1) where F > G, and false (0)
%   elsewhere.
%
% See also LT, GE, LE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

h = lt(g, f);

end