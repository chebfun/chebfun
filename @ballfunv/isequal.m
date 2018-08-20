function b=isequal(f, g)
% ISEQUAL Test the equality between two BALLFUNV
%   ISEQUAL(f, g) is the boolean f == g

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

b = iszero(f-g);
end
