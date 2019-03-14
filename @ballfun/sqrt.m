function g = sqrt(f)
%SQRT   Square root of a BALLFUN.
%   SQRT(F) is the square root of a BALLFUN F. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @sqrt ); 

end
