function out = rval(f)
%RVAL  Obtain the value of a ONEFUN at its right endpoint.
%   RVAL(F) returns the F(1).

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

out = feval(f, 1);

end