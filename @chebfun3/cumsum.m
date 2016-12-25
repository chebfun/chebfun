function f = cumsum(f, dim)
%CUMSUM   Indefinite integral of a CHEBFUN3.
%   F = CUMSUM(F) returns the indefinite integral of a CHEBFUN3 with 
%   respect to one variable and hence, returns a CHEBFUN3 object. The 
%   integration is done by default in the x-direction.
%
%   F = CUMSUM(F, DIM). If DIM = 1 integration is along the x-direction, if 
%   DIM = 2 integration is along the y-direction, and if DIM = 3 
%   integration is along the z-direction.
%
% See also CHEBFUN3/CUMSUM2, CHEBFUN3/CUMSUM3, CHEBFUN3/SUM, 
%   CHEBFUN3/SUM2 and CHEBFUN3/SUM3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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