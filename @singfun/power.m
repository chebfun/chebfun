function f = power(f, b)
% .^   SINGFUN power.
%   F.^G returns a SINGFUN F to the scalar power G, a scalar F to the SINGFUN
%   power G, or a SINGFUN F to the SINGFUN power G. F and or G may be complex.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

f.smoothPart = power(f.smoothPart, b);
f.exponents = b*f.exponents;

end