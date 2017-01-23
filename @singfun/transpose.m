function f = transpose(f) %#ok<*INUSD>
%TRANSPOSE   SINGFUN objects are not transposable.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SINGFUN:transpose:notPossible', ...
      'SINGFUN objects are not transposable.')

end
