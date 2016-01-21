function S = sum(disc)
%SUM   Definite integral operator for TRIGSPEC.
%   S = SUM(DISC) returns a definite integral operator for the TRIGSPEC
%   class.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Useful information:
dom = disc.domain;
N = disc.dimension;

% Rescale the differentiation matrix.
rescaleFactor = dom(2) - dom(1);
S = trigspec.diffmat(N) * (2*pi/rescaleFactor);

end