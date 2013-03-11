function f = uminus(f)
%UMINUS   Unary minus of a BNDFUN.
%   UMINUS(F) is the negative of F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Negate the onefun:
f.onefun = -f.onefun; 

end