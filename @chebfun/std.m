function out = std(f)
%STD    Standard deviation of a chebfun.
%   STD(F) is the standard deviation of the CHEBFUN F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% The standard deviation is the square root SQRT() of the variance VAR():
out = sqrt(var(f));

end