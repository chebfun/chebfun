function f = uminus(f)
%UMINUS   Unary minus of a DELTAFUN.
%   UMINUS(F) is the negative of F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Negate the FUNPART:
f.funPart = -f.funPart;

% Negate the impulses:
f.deltaMag = -f.deltaMag;

end
