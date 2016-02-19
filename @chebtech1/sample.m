function [points, values] = sample(f, n)
%SAMPLE   Sample a CHEBTECH1 at first-kind Chebyshev points.
%   [POINTS, VALUES] = SAMPLE(F, N) returns a vector POINTS of N first-kind
%   Chebyshev points and a vector VALUES of the values of F at those points.
%
%   [POINTS, VALUES] = SAMPLE(F) uses N = LENGTH(F).

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    if ( nargin < 2 )
        n = length(f);
    end

    points = chebtech1.chebpts(n);
    values = feval(f, points);
end
