function f = exp(f) 
%EXP  Exponential of a CHEBFUN3
%   EXP(F) returns the exponential of a CHEBFUN3. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    return 
end 

op = @(x,y,z) exp(feval(f, x, y, z)); % Resample.
f = chebfun3(op, f.domain, 'fiberDim', 3);          % Call constructor. 

end