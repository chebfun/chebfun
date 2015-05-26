function h = rem(x, y)
%REM   Remainder after division of two CHEBFUN objects.
%   REM(X, Y) returns X - n.*Y, where n = FIX(X./Y).
%
% See also FIX.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

n = fix(x./y);
h = x - n.*y;

end
