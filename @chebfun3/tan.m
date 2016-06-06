function f = tan(f)
%TAN   Tangent of a CHEBFUN3.
%   TAN(F) returns the tangent of a CHEBFUN3 object F.
%
% See also CHEBFUN3/TAND and CHEBFUN3/COMPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    return
end

op = @(x,y,z) tan(feval(f, x, y, z));   % Resample
f = chebfun3(op, f.domain);             % Call constructor.

end