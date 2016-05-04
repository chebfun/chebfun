function F = uminus( F )
%UMINUS   Unary minus for a SEPARABLEAPPROX. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

F.pivotValues = -F.pivotValues;

end
