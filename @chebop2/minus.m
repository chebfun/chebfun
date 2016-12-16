function N = minus(N1, N2)
%MINUS   Minus for CHEBOP2 objects.
%
% N = MINUS(N1, N2) is the same as N = N1 - N2. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

N = plus(N1, uminus(N2));

end
