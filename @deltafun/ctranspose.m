function f = ctranspose(f) %#ok<*INUSD>
%CTRANSPOSE   Conjugate transpose of a DELTAFUN object.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = conj(f);
f.isTransposed = ~f.isTransposed;
%[TODO]: What should we do with the funPart of f?
end
