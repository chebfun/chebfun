function val = lval(fun)
%LVAL  Obtain the value of a FUN at its left endpoint.
%   LVAL(F) returns the value of F at F.domain(1).

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call LVAL() of the FUNPART:
val = lval(fun.funPart);

end