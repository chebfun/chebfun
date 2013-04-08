function transpose(f) %#ok<*INUSD>
%TRANSPOSE   CHEBTECH objects are not transposable.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

    error('CHEBFUN:CHEBTECH:transpose:notpossible', ...
        'CHEBTECH objects are not transposable.')
    
end
