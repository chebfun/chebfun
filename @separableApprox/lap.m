function L = lap(f)
%LAP   Laplacian of a SEPARABLEAPPROX.
%   L = LAP(F) returns a SEPARABLEAPPROX representing the Laplacian of F. 
%
%   This is shorthand for LAPLACIAN(F).
%
% See also SEPARABLEAPPROX/LAPLACIAN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call SEPARABLEAPPROX/LAPLACIAN:
L = laplacian(f);

end
