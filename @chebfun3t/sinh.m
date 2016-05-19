function f = sinh(f)
%SINH   Hyperbolic sine of a CHEBFUN3T object.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) )
    return
end 

op = @(x,y,z) sinh(feval(f, x, y, z));     % Resample.
f = chebfun3t(op, f.domain);               % Call constructor. 

end