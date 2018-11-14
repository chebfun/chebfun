function g = real(f)
%REAL Real part of a BALLFUN function
%   REAL(F) is the real part of the BALLFUN F.
%
% See also IMAG. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @real ); 

end
