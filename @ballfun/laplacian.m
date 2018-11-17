function g = laplacian(f)
%LAPLACIAN   Laplacian of a BALLFUN.
%   L = LAPLACIAN(F) returns a BALLFUN representing the Laplacian of F.
%
% See also CHEBFUN2/LAPLACIAN.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = diff(f,1,2) + diff(f,2,2) + diff(f,3,2);

end
