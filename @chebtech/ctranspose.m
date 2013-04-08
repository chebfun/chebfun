function ctranspose(f) %#ok<*INUSD>
%CTRANSPOSE  CHEBTECH objects are not transposable.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

    error('CHEBFUN:CHEBTECH:ctranspose:notpossible', ...
        'CHEBTECH objects are not transposable.')
    
end
