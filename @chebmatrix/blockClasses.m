function t = blockClasses(L)
%BLOCKCLASSES Class of each block in the chebmatrix.
%   Returns a cell arrray. 

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

t = cellfun(@class, L.blocks, 'uniform', false);

end
