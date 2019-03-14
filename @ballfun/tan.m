function g = tan(f)
%TAN   Tangent of a BALLFUN.
%   TAN(F) computes the tangent of the BALLFUN F.
%
% See also TANH, SIN, COS.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @tan ); 
end
