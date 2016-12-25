function f = cumsum(f, dim)
%CUMSUM   Indefinite integral of a SEPARABLEAPPROX.
%   F = CUMSUM(F) returns the indefinite integral of a SEPARABLEAPPROX with respect to
%   one variable and hence, returns a chebfun. The integration is done by
%   default in the y-direction.
%
%   F = CUMSUM(F, DIM). If DIM = 2 integration is along the x-direction, if DIM
%   = 1 integration is along the y-direction.
%
% See also CUMSUM2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty( f ) )  
    f = [];
    return
end

if ( nargin == 1 )
    % Default to integration along the y-direction.
    dim = 1;
end

if ( dim == 1 )
    % CUMSUM along the y-variable.
    f.cols = cumsum( f.cols );
elseif ( dim == 2 )
    % CUMSUM along the rows.
    f.rows = cumsum( f.rows );
else
    error('CHEBFUN:SEPARABLEAPPROX:cumsum:dim', ...
        'Integration direction must be x or y.');
end

end
