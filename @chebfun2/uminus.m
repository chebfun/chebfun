function F = uminus( F )
%UMINUS   Unary minus for a CHEBFUN2. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

F.pivotValues = -F.pivotValues;

end
