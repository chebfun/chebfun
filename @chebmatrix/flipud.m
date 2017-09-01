function F = flipud(F)
%FLIPUD   Flip the rows of a chebmatrix.

%  Copyright 2017 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

F.blocks = flipud(F.blocks);

end
