function vals = feval(F, x, y, z)
%FEVAL   Pointwise evaluation of a CHEBFUN3V object.
%   F(X, Y, Z) returns value of F at the point (X,Y,Z).
%
% See also CHEBFUN3V/SUBSREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) )
    vals = []; 
    return
end

nF = F.nComponents;
vals = zeros(nF, length(x));

% Evaluate each component:
for jj = 1:nF
   vals(jj, :) = feval(F.components{jj}, x, y, z);  
end

end