function f = uminus(f)
%UMINUS   Unary minus of a CLASSICFUN.
%   UMINUS(F) is the negative of F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Negate the ONEFUN:
f.onefun = -f.onefun; 

end
