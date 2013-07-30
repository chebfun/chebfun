function g = sinc(f)
%SINC   Sinc function of a chebfun.
%
% See also SIN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

g = comp(f, @(x) sin(x)./x);

end