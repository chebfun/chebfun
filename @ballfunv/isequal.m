function b=isequal(f, g)
%ISEQUAL   Equality test for two BALLFUNVs.
%   ISEQUAL(f, g) returns logical 1 (TRUE) if the BALLFUNV objects F and G
%   contain identical components, and logical 0 (FALSE) otherwise.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

b = iszero(f-g);
end
