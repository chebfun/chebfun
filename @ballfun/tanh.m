function g = tanh(f)
%TANH   Hyperbolic tangent of a BALLFUN.
%   TANH(F) computes the hyperbolic tangent of the BALLFUN F.
%
% See also TAN, SINH, COSH.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @tanh ); 
end
