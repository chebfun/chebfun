function g = sinh(f)
%SINH   Hyperbolic sine of a BALLFUN.
%   SINH(F) computes the hyperbolic sine of the BALLFUN F.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @sinh ); 
end
