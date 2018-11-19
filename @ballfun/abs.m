function g = abs( f )
%ABS Absolute value of a BALLFUN.
%   ABS(F) is the absolute value of the BALLFUN F. This function does not work
%   if the function passes through or becomes numerically close to zero.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @abs ); 
end
