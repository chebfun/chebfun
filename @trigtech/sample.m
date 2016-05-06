function [values, points] = sample(f, n)
%SAMPLE   Sample a TRIGTECH at equally-spaced points.
%   VALUES = SAMPLE(F, N) returns a vector VALUES of the values of F at N
%   equally-spaced points.
%
%   [VALUES, POINTS] = SAMPLE(F, N) returns also the vector POINTS at which the
%   values were computed.
%
%   [...] = SAMPLE(F) uses N = LENGTH(F).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    if ( nargin < 2 )
        n = length(f);
    end

    if ( n == length(f) )
        values = f.values;
    else
        values = trigtech.coeffs2vals(trigtech.alias(f.coeffs, n));
    end

    if ( nargout > 1 )
        points = trigtech.trigpts(n);
    end
end
