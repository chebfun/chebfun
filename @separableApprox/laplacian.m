function L = laplacian(f)
%LAPLACIAN   Laplacian of a SEPARABLEAPPROX.
%   L = LAPLACIAN(F) returns a SEPARABLEAPPROX representing the Laplacian of F.
%
% See also SEPARABLEAPPROX/LAP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% laplacian(f) = f_xx + f_yy: 
L = diff(f, 2, 2) + diff(f, 2, 1); 

end
