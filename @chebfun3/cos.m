function f = cos(f)
%COS   Cosine of a CHEBFUN3.
%   COS(F) returns the cosine of F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(f) )
    return
end 

op = @(x,y, z) cos(feval(f, x, y, z));  % Resample. 
f = chebfun3(op, f.domain, 'fiberDim', 3);          % Call constructor. 

end