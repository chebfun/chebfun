function r = range(f, dim)
%RANGE   Range of CHEBFUN.
%   R = RANGE(F) returns the range R = MAX(F) - MIN(F) of the CHEBFUN F.
%
%   R = RANGE(F, DIM) operates along the dimension DIM of the quasimatrix F. If
%   DIM represents the continuous variable, then R is a vector. If DIM
%   represents the discrete dimension, then R is a CHEBFUN. The default for
%   DIM is 1, unless F has a singleton dimension, in which case DIM is the
%   continuous variable.
%
% See also CHEBFUN/MINANDMAX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    r = [];
    return
end

if ( nargin < 2 )
    if ( numColumns(f) > 1 )
        % Dim 1 by default.
        dim = 1;
    else
        % Along continuous dimension for scalar f.
        dim = 1 + f(1).isTransposed;
    end
end

if ( dim >= 3 )
    % Consistent with MATLAB?
    r = 0*f;
    return
end

% Compute the MIN and MAX:
r = minandmax(f, [], dim);

if ( isnumeric(r) )
    % Range along continuous dimension.
    r = diff(r);

    % Handle row inputs:
    if ( f(1).isTransposed )
        r = r.';
    end
else
    % Range across discrete dimension.
    if ( f(1).isTransposed )
        r = diff(extractColumns(r, 1:2), 1, 1);
    else
        r = diff(extractColumns(r, 1:2), 1, 2);
    end
end

end
