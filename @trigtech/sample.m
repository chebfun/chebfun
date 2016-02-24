function [values, points] = sample(f, n)
%SAMPLE   Sample a TRIGTECH at equally-spaced points.
%   [VALUES, POINTS] = SAMPLE(F, N) returns a vector POINTS of N equally-spaced
%   points and a vector VALUES of the values of F at those points.
%
%   [VALUES, POINTS] = SAMPLE(F) uses N = LENGTH(F).

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    if ( nargin < 2 )
        n = length(f);
    end

    points = trigtech.trigpts(n);
    if ( n == length(f) )
        values = f.values;
    else
        values = trigtech.coeffs2vals(trigtech.alias(f.coeffs, n));
    end
end
