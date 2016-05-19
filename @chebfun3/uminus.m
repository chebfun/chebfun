function F = uminus(F)
%UMINUS   Unary minus for CHEBFUN3 objects.
%
% See also CHEBFUN3/UPLUS and CHEBFUN3/MINUS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F.core = -F.core;

end