function h = minus(f, g)
%-   BALLFUNV minus.
%   F - G subtracts BALLFUNVs F and G, or a scalar from a BALLFUNV if either F or
%   G is a scalar.
%
%   H = MINUS(F, G) is called for the syntax 'F - G'.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

h = f + (-g);
end
