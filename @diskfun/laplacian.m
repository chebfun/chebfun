function f = laplacian( f ) 
% LAPLACIAN     Scalar laplacian of a diskfun.
%   L = LAPCIAN(F) returns a DISKFUN representing the Laplacian of F. 
%
% See also DISKFUN/DIFF, DISKFUN/GRADIENT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = diff(f, 1, 2) + diff( f, 2, 2); 

end 