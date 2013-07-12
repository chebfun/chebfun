function ctranspose(f) %#ok<*INUSD>
%CTRANSPOSE  SINGFUN objects are not transposable, so this method will throw an
%   error.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

error('CHEBFUN:SINGFUN:ctranspose:notpossible', ...
    'SINGFUN objects are not transposable.')
    
end
