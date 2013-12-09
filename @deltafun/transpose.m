function f = transpose(f) %#ok<*INUSD>
%TRANSPOSE   Transpose a DELTAFUN object.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f.isTransposed = ~f.isTransposed;
%[TODO]: What should we do with the funPart of f?
end
