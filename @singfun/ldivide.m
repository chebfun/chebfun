function f = ldivide(f, g)
%.\   Left divide for a SINGFUN.
%   F .\ G is equivalent to G ./ F.
%
% See also RDIVIDE, TIMES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call RDIVIDE():
f = rdivide(g, f);

end
