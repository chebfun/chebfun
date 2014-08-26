function D = diff(disc, m)
%DIFF    Differentiation operator for FOURCOLLOC discretization.
%   D = DIFF(DISC) gives the matrix such that if v=D*u, then v=u', where u
%   is a FOURCOLLOC representation of a trigonometric polynomial.
%
%   DIFF(DISC, M) for positive integer M returns D^M.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

% Store information about domain and dimensions.
d = disc.domain;
n = disc.dimension;

if ( m == 0 )
    % Trivial case.
    D = eye(sum(n));
else
    len = d(2) - d(1);
    % Rescale the differentiation matrix.
    D = fourcolloc.diffmat(n, m) * (2*pi/len)^m;
end