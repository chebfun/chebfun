function f = ldivide(f, g)
%.\   Left array divide for a DELTAFUN.
%   F .\ G is equivalent to G ./ F.
%
% See also MLDIVIDE, RDIVIDE, TIMES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = rdivide(g, f);

end
