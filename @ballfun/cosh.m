function g = cosh(f)
%COSH   Hyperbolic cosine of a BALLFUN.
%   COSH(F) computes the hyperbolic cosine of F.
%
% See also SINH, COS, SIN

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @cosh ); 
end
