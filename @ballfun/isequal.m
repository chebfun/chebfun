function b = isequal(f, g)
%ISEQUAL Test the equality between two BALLFUNs
%   ISEQUAL(F, G) returns logical 1 (TRUE) if the BALLFUN objects F and G
%   are identical, and logical 0 (FALSE) otherwise.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Test if f = g
b = iszero( f-g );
end
