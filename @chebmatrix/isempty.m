function t = isempty(L)
%ISEMPTY(A)   True if there are no blocks in the CHEBMATRIX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

t = isempty(L.blocks);

end
