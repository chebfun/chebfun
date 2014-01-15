function d = dirac(f)
%DIRAC   Dirac delta function.
%
% See also HEAVISIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

f = 0*f;
d = deltafun( [], 1, 0 );
