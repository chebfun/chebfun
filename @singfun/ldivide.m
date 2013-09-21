function f = ldivide(f, g)
%.\   Left array divide for a SINGFUN.
%   F .\ G is equivalent to G ./ F.
%
% See also RDIVIDE, TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = rdivide(g, f);

end
