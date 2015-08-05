function f = ctranspose(f) %#ok<*INUSD>
%CTRANSPOSE   SINGFUN objects are not transposable.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SINGFUN:ctranspose:notPossible', ...
      'SINGFUN objects are not transposable.')
    
end
