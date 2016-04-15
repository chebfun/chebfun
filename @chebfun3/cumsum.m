function f = cumsum(f, dim)
%CUMSUM   Indefinite integral of a CHEBFUN3.
%   F = CUMSUM(F) returns the indefinite integral of a CHEBFUN3 with 
%   respect to one variable and hence, returns a chebfun3 object. The 
%   integration is done by default in the x-direction.
%
%   F = CUMSUM(F, DIM). If DIM = 1 integration is along the x-direction, if 
%   DIM = 2 integration is along the y-direction, and if DIM = 3 
%   integration is along the z-direction.
%
%   See also chebfun3/cumsum2, chebfun3/cumsum3, chebfun3/sum, 
%   chebfun3/sum and chebfun3/sum3.

% Check for empty:
if ( isempty(f) )
    f = [];
    return
end

if ( nargin == 1 )
    % Default to integration along the x-direction.
    dim = 1;
end

if ( dim == 1 )
    % CUMSUM along the 1st variable.
    f.cols = cumsum(f.cols);
elseif ( dim == 2 )
    % CUMSUM along the 2nd variable.
    f.rows = cumsum(f.rows);
elseif ( dim == 3 )    
    % CUMSUM along the 3rd variable.
    f.tubes = cumsum(f.tubes);
else
    error('CHEBFUN:CHEBFUN3:cumsum:dim', ...
        'Integration direction must be x, y, or z.');
end

end