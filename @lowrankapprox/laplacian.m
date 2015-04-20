function L = laplacian(f)
%LAPLACIAN   Laplacian of a CHEBFUN2.
%   L = LAPLACIAN(F) returns a CHEBFUN2 representing the Laplacian of F.
%
% See also LAP.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% laplacian(f) = f_xx + f_yy: 
L = diff(f, 2, 2) + diff(f, 2, 1); 

end
