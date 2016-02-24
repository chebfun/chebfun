function [values, points] = sample(f, varargin{:})
%SAMPLE   Sample a SINGFUN on an "appropriate" grid.
%   [VALUES, POINTS] = SAMPLE(F, N) returns a vector POINTS of N
%   "appropriately-chosen" points and a vector VALUES of the values of F at
%   those points.  What "appropriately-chosen" means depends on the type of
%   F.SMOOTHPART.
%
%   [VALUES, POINTS] = SAMPLE(F) uses N = LENGTH(F).

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [values, points] = sample(f.smoothPart);

    % Multiply now with the singular factors:
    if ( f.exponents(1) )
        values = values.*(1 + x).^f.exponents(1);
    end

    if ( f.exponents(2) )
        values = values.*(1 - x).^f.exponents(2);
    end
end
