function L = del2( f )
%DEL2   Scaled Laplacian of a CHEBFUN2.
%   L = del2(F) is a numerical approximation of 
%       del^2 F/4 = (d^2F/dx^2 + d^2F/dy^2)/4.
%
% See also LAPLACIAN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

L = ( diff(f, 2 , 2) + diff(f, 2, 1) ) / 4; 

end