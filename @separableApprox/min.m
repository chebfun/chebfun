function g = min( f, g, dim )
%MIN   Minimum value of a SEPARABLEAPPROX in one direction.
%   MIN(f) returns a chebfun representing the minimum of the SEPARABLEAPPROX along the
%   y direction, i.e, MIN(f) = @(x) max( f ( x, : ) )
%
%   MIN(f, [], dim) returns a CHEBFUN representing the minimum of f along the
%   DIM direction. If DIM = 1 is along the y-direction and DIM = 2 is along the
%   x-direction.
%
%   WARNING: This function is not always accurate to the expected precision.
% 
%   For the global minimum use MIN2.
%
% See also MAX, MAX2, MIN2, MINANDMAX2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    error('CHEBFUN:SEPARABLEAPPROX:min:input', 'SEPARABLEAPPROX is empty');
end

% Default to max of one chebfun2:
if ( nargin < 2 )
    g = []; 
end

% Default to maximum along the y direction: 
if ( nargin < 3 )
    dim = 1;
end

% Do not allow min(F, G).
if ( nargin > 1 && ~isempty(g) )   
    error('CHEBFUN:SEPARABLEAPPROX:min:twoSeparableApproxInputs', ...
        'Unable to minimise two SEPARABLEAPPROX objects.');
end

% min(f) = - max ( -f )  
g = -max( -f, [], dim );

end
