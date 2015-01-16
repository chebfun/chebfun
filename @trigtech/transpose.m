function transpose(f)
%TRANSPOSE   TRIGTECH objects are not transposable, so this method will throw an
%error.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:TRIGTECH:transpose:notpossible', ...
    'TRIGTECH objects are not transposable.')
    
end
