function f = cumsum(f,dim)
%CUMSUM Indefinite integral of a chebfun2.
%
% F = CUMSUM(F) returns the indefinite integral of a chebfun2 with respect
% to one variable and hence, returns a chebfun. The integration is done
% by default in the y-direction.
%
% F = CUMSUM(F, DIM). If DIM = 2 integration is along the x-direction, if
% DIM = 1 integration is along the y-direction.
%
% See also CUMSUM2.

if ( nargin == 1 ) % default to integration along the y-direction.
    dim = 1;
end

if ( isempty( f ) )  % check for empty chebfun2.
    f = [];
    return
end

if ( dim == 1 )
    % cumsum along the y-variable.
    f.cols = cumsum( f.cols );
elseif ( dim == 2 )
    % cumsum along the rows.
    f.rows = cumsum( f.rows );
else
    error('CHEBFUN2:CUMSUM:DIM','Integration direction must be x or y');
end

end