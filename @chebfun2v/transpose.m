function F = transpose( F )
% .' transpose of a CHEBFUN2V

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

F.isTransposed = ~F.isTransposed;

end
