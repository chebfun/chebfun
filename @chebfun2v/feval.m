function vals = feval(F, x, y)
%FEVAL pointwise evaluate a chebfun2v.
%   F(X,Y) returns the evaluation of F at the coordinate (X,Y).
%
% See also SUBSREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( F ) ) 
    vals = []; 
    return
end

nF = F.nComponents; 
vals = zeros(nF, length(x)); 

% Evaluate each component:
for jj = 1:nF
   vals(jj, :) = feval(F.components{jj}, x, y);  
end


end