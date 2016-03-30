function f = minus(f, g)
%-   SPHEREFUNV minus.
%   F - G subtracts the SPHEREFUNV F from G componentwise.
%
%   MINUS(F, G) is called for the syntax 'F - G'.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = plus(f, uminus(g));

end
