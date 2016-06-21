function f = tand(f)
%TAND  Tangent of a CHEBFUN3 in degrees.
%   TAND(F) returns the tangent of a CHEBFUN3 object F in degrees.
%
% See also CHEBFUN3/TAN and CHEBFUN3/COMPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) )
    return
end

op = @(x,y,z) tand(feval(f, x, y, z));      % Resample.
f = chebfun3(op, f.domain);                 % Call constructor.

end