function f = sin(f)
%SIN   Sine of a CHEBFUN3.
%
%   SIN(F) returns the sine of F.
%
%   See also CHEBFUN3/COMPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(f) )
    return
end 

op = @(x,y, z) sin(feval(f, x, y, z));  % Resample. 
f = chebfun3(op, f.domain, 'fiberDim', 3);           % Call constructor. 

end