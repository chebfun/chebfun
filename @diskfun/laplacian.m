function f = laplacian( f ) 
% LAPLACIAN     Scalar laplacian of a diskfun.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = diff(f, 1, 2) + diff( f, 2, 2); 

end 