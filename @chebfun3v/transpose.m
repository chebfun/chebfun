function F = transpose(F)
%.'   Transpose of a CHEBFUN3V object.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F.isTransposed = ~F.isTransposed;

end