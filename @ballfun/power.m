function g = power(f,n)
%.^   BALLFUN power.
%   F.^G returns a BALLFUN F to the scalar power G.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = ballfun( n );
g = compose(f, @power, n); 
end
