function vals = feval(F, x, y, z)
%FEVAL pointwise evaluate a CHEBFUN3V.
%
%   F(X,Y) returns the evaluation of F at the coordinate (X,Y,Z).
%
%   See also CHEBFUN3V/SUBSREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) ) 
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