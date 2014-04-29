function f = transpose(f) %#ok<*INUSD>
%TRANSPOSE   SINGFUN objects are not transposable.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SINGFUN:transpose:notpossible', ...
      'SINGFUN objects are not transposable.')

end
