function out = cov(f)
%COV   Covariance of a chebfun.
%   COV(F) is the covariance of the chebfun F. This is the same as the variance.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% The covariance f a chebfun is the same as the variance:
out = var(f);

% [TODO]: Not for array-valued chenfun objects!

end

