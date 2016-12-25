function f = sin(f)
%SIN   Sine of a CHEBFUN3.
%   SIN(F) returns the sine of a CHEBFUN3 object F.
%
% See also CHEBFUN3/COMPOSE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(f) )
    return
end 

op = @(x,y,z) sin(feval(f, x, y, z));   % Resample.
f = chebfun3(op, f.domain);             % Call constructor.

end