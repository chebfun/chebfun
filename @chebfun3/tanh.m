function f = tanh(f)
%TANH   Hyperbolic tangent of a CHEBFUN3.
%   TANH(F) returns the hyperbolic tangent of a CHEBFUN3 object F.
%
% See also CHEBFUN3/COMPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    return
end

op = @(x,y,z) tanh(feval(f, x, y, z));  % Resample.
f = chebfun3(op, f.domain);             % Call constructor.

end