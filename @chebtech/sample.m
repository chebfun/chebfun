function [values, points] = sample(f, n)
%SAMPLE   Sample a CHEBTECH at appropriate Chebyshev points.
%   [VALUES, POINTS] = SAMPLE(F, N) returns a vector POINTS of N first- or
%   second-kind Chebyshev points (according to the type of F) and a vector
%   VALUES of the values of F at those points.
%
%   [VALUES, POINTS] = SAMPLE(F) uses N = LENGTH(F).

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    if ( nargin < 2 )
        n = length(f);
    end

    points = f.chebpts(n);
    values = f.coeffs2vals(f.alias(f.coeffs, n));
end

