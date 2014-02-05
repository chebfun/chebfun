function g = min( f, g, dim )
%MIN  Minimum value of a CHEBFUN2.
%   MIN(F), returns a CHEBFUN representing the minimum of the CHEBFUN2 along the
%   y direction, i.e., MIN(F) = @(x) min( F ( x, : ) )
% 
%   MIN(F,[],DIM) returns a CHEBFUN representing the minimum of f along the DIM
%   direction.
% 
% See also MAX. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    error('CHEBFUN2:MIN:INPUT','Chebfun2 is empty');
end

% Do not allow min(F, G).
if ( nargin > 1 && ~empty(g) )   
    error('CHEBFUN2:MIN', 'Unable to minimise two CHEBFUN2 objects.');
end

% Default to maximum along the y direction:
if ( nargin < 3 )   
    dim = 1;
end

% min(f) = - max ( -f )  
g = -max( -f, [], dim );

end