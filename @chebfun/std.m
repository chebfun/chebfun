function out = std(f)
%STD   Standard deviation of a CHEBFUN.
%   STD(F) is the standard deviation of the CHEBFUN F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% The standard deviation is the square root SQRT() of the variance VAR():
out = sqrt(var(f));

end
