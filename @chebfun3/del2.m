function L = del2(f)
%DEL2   Scaled Laplacian of a CHEBFUN3.
%   L = del2(F) is a Chebfun3 representation of
%       del^2 F/6 = (d^2F/dx^2 + d^2F/dy^2 + d^2F/dz^2)/6.
%
% See also LAPLACIAN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

L = laplacian(f)/6;

end
