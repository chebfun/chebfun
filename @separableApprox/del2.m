function L = del2( f )
%DEL2   Scaled Laplacian of a SEPARABLEAPPROX.
%   L = del2(F) is a numerical approximation of 
%       del^2 F/4 = (d^2F/dx^2 + d^2F/dy^2)/4.
%
% See also LAPLACIAN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

L = laplacian( f )/4;

end
