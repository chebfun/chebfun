function D = diff(disc, m)
%DIFF    Differentiation operator for TRIGCOLLOC discretization.
%   D = DIFF(DISC) gives the matrix such that if v=D*u, then v=u', where u
%   is a TRIGCOLLOC representation of a trigonometric polynomial.
%
%   DIFF(DISC, M) for positive integer M returns D^M.

%  Copyright 2016 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

% Store information about domain and dimensions.
dom = disc.domain;
n = disc.dimension;

if ( m == 0 )
    % Trivial case.
    D = eye(sum(n));
else
    rescaleFactor = dom(2) - dom(1);
    % Rescale the differentiation matrix.
    D = trigcolloc.diffmat(n, m) * (2/rescaleFactor)^m;
end

end
