function F = uminus(F)
%UMINUS   Unary minus for CHEBFUN3 objects.
%   -F negates the CHEBFUN3 object F.
%
%   G = UMINUS(F) is called for the syntax '-F'.
%
% See also CHEBFUN3/UPLUS and CHEBFUN3/MINUS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F.core = -F.core;

end