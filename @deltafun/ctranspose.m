function f = ctranspose(f) %#ok<*INUSD>
%CTRANSPOSE   DELTAFUN objects are not transposable.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:DELTAFUN:ctranspose:notPossible', ...
    'DELTAFUN objects are not transposable.')

end
