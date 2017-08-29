function F = fliplr(F)
%FLIPLR   Flip the columns of a CHEBMATRIX.

%  Copyright 2017 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

F.blocks = fliplr(F.blocks);

end
