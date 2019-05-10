function b=isequal(f, g)
%ISEQUAL Equality test for BALLFUNV.  
%   ISEQUAL(F, G) returns 1 then F and G are the same
%   BALLFUNV and 0 otherwise.
%
%   See also ISZERO.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

b = iszero(f-g);
end
