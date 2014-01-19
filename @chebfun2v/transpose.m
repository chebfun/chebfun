function F = transpose(F)
% .' transpose of a chebfun2v

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

F.isTransposed = ~F.isTransposed;

end