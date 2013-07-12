function f = uminus(f)
%UMINUS Negate a singfun.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

f.smoothPart = -f.smoothPart;

end