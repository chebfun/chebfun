function transpose(f) %#ok<*INUSD>
%TRANSPOSE   CLASSICFUN objects are not transposable.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:CLASSICFUN:transpose:notpossible', ...
    'CLASSICFUN objects are not transposable.')

end
