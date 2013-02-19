function F = mldivide(c, f)
%\	Left scalar divide of a FUNCHEB1.
%   C \ F divides the FUNCHEB1 F by a scalar C.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO] Add support for vectorised FUNCHEB1 objects.

% Simply call RDIVIDE (since C is a scalar):
F = mrdivide(f, c);

end
