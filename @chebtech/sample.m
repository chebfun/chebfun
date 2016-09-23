function [values, points] = sample(f, n)
%SAMPLE   Sample a CHEBTECH at appropriate Chebyshev points.
%   VALUES = SAMPLE(F, N) returns a vector VALUES of the values of F at N
%   first- or second-kind Chebyshev points (according to the type of F).
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

    values = f.coeffs2vals(f.alias(f.coeffs, n));

    if ( nargout > 1 )
        points = f.chebpts(n);
    end
end
