function g = fracInt(f, mu)
%FRACINT  Fractional integral of a SINGFUN. 
%   FRACINT(F, MU) gives the order MU fractional integral of a SINGFUN object F.
%
%   Currently this only supports the situation where F is smooth (i.e., it has
%   trivial exponents).

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( any( f.exponents ) )
    error('CHEBFUN:SINGFUN:fracInt:exponents', ...
        ['FRACINT currently this only supports the situation where g is smooth',
         ' (i.e., it has trivial exponents)']);
end

g = f;
g.smoothPart = fracInt(f.smoothPart, mu);
g.exponents = [mu, 0];

end