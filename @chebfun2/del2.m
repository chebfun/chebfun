function L = del2( f )
%DEL2 Scaled Laplacian of a chebfun2.
% 
% L = del2(f) when f is a chebfun2 is a numerical approximation of 
% 0.25*del^2 f = (d^2f/dx^2 + d^2f/dy^2)/4, where L is a chebfun2. 
%
% See also LAPLACIAN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

L = ( diff(f, 2 , 2) + diff(f, 2, 1) ) / 4; 

end