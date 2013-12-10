function f = uminus(f)
%UMINUS   Unary minus of a DELTAFUN.
%   UMINUS(F) is the negative of F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

f.funPart = -f.funPart;
f.impulses = -f.impulses;

end