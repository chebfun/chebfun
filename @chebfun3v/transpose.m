function F = transpose( F )
% .' transpose of a CHEBFUN3V

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F.isTransposed = ~F.isTransposed;

end