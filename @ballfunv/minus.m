function h = minus(f, g)
%-   BALLFUNV minus.
%   F - G subtracts the BALLFUNV F from G componentwise.
%
%   H = MINUS(F, G) is called for the syntax 'F - G'.
% 
% See also PLUS.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

h = f + (-g);
end
