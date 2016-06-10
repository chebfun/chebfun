function L = lap( f )
%LAP   Laplacian of a CHEBFUN3V.
%
%   L = LAP(F) returns a CHEBFUN3V representing the Laplacian of F. 
%
%   This is shorthand for LAPLACIAN( F )
%
%   See also LAPLACIAN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call Laplacian: 
L = laplacian( f );

end