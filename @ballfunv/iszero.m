function b = iszero(f)
% ISZERO Test the nullity of a BALLFUNV up to machine precision
%   ISZERO(f) is the boolean |f1| < 1e-10 && |f2| < 1e-10 && |f3| < 1e-10

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.comp;
b = (iszero(F{1}) && iszero(F{2}) && iszero(F{3}));
end
