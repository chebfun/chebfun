function transpose(f) %#ok<*INUSD>
%TRANSPOSE   CHEBTECH objects are not transposable, so this method will throw an
%   error.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:CHEBTECH:transpose:notPossible', ...
    'CHEBTECH objects are not transposable.')
    
end
