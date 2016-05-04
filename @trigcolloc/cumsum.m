function C = cumsum(disc, m)
%CUMSUM   Indefinite integration operator for TRIGCOLLOC discretization.
%   C = CUMSUM(DISC) gives the matrix such that if v=C*u, then u=v' and v=0 at
%   the left endpoint, as accurately as possible in trigonometric polynomial
%   discretization.
%
%   CUMSUM(DISC, M) for positive integer M returns C^M.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Store information about domain and dimensions.
dom = disc.domain;
n = disc.dimension;

if ( m == 0 )
    % Trivial case
    C = eye(sum(n));
else
    rescaleFactor = dom(2) - dom(1);
    % Rescale.
    C = trigcolloc.cumsummat(n) * (rescaleFactor/(2*pi)); 
    C = C^m;
end

end
