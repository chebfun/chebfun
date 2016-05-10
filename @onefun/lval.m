function out = lval(f)
%LVAL  Obtain the value of a ONEFUN at its left endpoint.
%   LVAL(F) returns F(-1).

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

out = feval(f, -1);

end