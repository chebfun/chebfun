function f = sin(f)
%SIN   Sine of a CHEBFUN3T.
%   SIN(F) returns the sine of a CHEBFUN3T object F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(f) )
    return
end 

op = @(x,y,z) sin(feval(f, x, y, z));   % Resample. 
f = chebfun3t(op, f.domain);            % Call constructor. 

end