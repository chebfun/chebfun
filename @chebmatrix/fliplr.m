function F = fliplr(F)
%FLIPLR   Flip the columns of a CHEBMATRIX.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

F.blocks = fliplr(F.blocks);

end
