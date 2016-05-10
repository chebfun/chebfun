function C = cumsum(disc, m)
%CUMSUM   Indefinite integration operator for CHEBCOLLOC discretization.
%   C = CUMSUM(DISC) gives the matrix such that if v=C*u, then u=v' and v=0
%   at the left endpoint, as accurately as possible in Chebyshev polynomial
%   discretization.
%
%   CUMSUM(DISC, M) for positive integer M returns C^M.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Store information about domain and dimensions.
d = disc.domain;
n = disc.dimension;

if ( m == 0 )
    % Trivial case:
    C = eye(sum(n));
else
    numIntervals = disc.numIntervals;
    
    % Find the diagonal blocks.
    blocks = cell(numIntervals);
    for k = 1:numIntervals
        len = d(k+1) - d(k);
        blocks{k} = disc.cumsummat(n(k)) * (len/2);  % Scaled cumsummats.
    end
    
    % Assemble!
    C = blkdiag(blocks{:});
    
    % Each subinterval also contributes to the integrals in all the
    % subintervals to its right, creating a triangular structure.
    offset = 0;
    for k = 1:numIntervals
        % Grab the weights for the integral using all of this
        % subinterval.
        row = offset + n(k);
        cols = offset + (1:n(k));
        last = C(row, cols);
        % Copy it to add to the ones that follow.
        offset = row;
        C(offset + 1:end, cols) = repmat(last,[sum(n) - offset, 1]);
    end
    C = C^m;
end

end

