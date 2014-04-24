function I = eye(disc)
%EYE   Identity operator for COLLOC discretization.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = disc.dimension;
I = eye(sum(n));

end
