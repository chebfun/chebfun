function L = lap( F )
%LAP   Vector Laplacian of a CHEBFUN2V.
%   LAP( F ) returns a CHEBFUN2 representing the Laplacian of F. 
%
%   This is shorthand for LAPLACIAN( F )
%
% See also LAPLACIAN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

L = laplacian( F );

end