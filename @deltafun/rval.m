function val = rval(fun)
%RVAL  Obtain the value of a FUN at its right endpoint.
%   RVAL(F) returns the value of F at F.domain(end).

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call RLVAL() of the ONEFUN:
val = rval(fun.funPart);

end