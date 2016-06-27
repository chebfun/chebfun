function f = cosh(f)
%COSH   Hyperbolic cosine of a CHEBFUN3.
%   COSH(F) returns the hyperbolic cosine of a CHEBFUN3 object F.
%
% See also CHEBFUN3/COMPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(f) ) 
    return 
end 

op = @(x,y,z) cosh(feval(f, x, y, z));  % Resample.
f = chebfun3(op, f.domain);             % Call constructor.

end