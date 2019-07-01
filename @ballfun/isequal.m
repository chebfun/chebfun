function b = isequal(f, g)
%ISEQUAL Equality test for BALLFUN.  
%   ISEQUAL(F, G) returns 1 if F =G, 
%   returns 0 otherwise 
%
%   See also ISZERO.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Test if f = g
b = iszero( f - g );
end
