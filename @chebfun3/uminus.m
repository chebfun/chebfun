function F = uminus(F)
%UMINUS   Unary minus for a CHEBFUN3. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F.core = -F.core;

end