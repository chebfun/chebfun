function f = cosh(f)
%COSH   Hyperbolic cosine of a CHEBFUN3T object.
%
%   COSH(F) returns the hyperbolic cosine of F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(f) )
    return 
end 

op = @(x,y,z) cosh(feval(f, x, y, z)) ;  % Resample.
f = chebfun3t(op, f.domain);             % Call constructor.

end