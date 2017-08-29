function L = lap(f)
%LAP   Vector Laplacian of a CHEBFUN3V object.
%   L = LAP(F) returns a CHEBFUN3V representing the vector Laplacian of F.
%
%   This is shorthand for LAPLACIAN(F).
%
% See also CHEBFUN3V/LAPLACIAN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call Laplacian: 
L = laplacian(f);

end