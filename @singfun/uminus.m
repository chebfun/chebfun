function f = uminus(f)
%UMINUS Negate a SINGFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

f.smoothPart = -f.smoothPart;

end