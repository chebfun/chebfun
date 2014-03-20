function transpose(f)
%TRANSPOSE   FOURIETECH objects are not transposable, so this method will throw an
%   error.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

error('FOURIETECH:FOURIETECH:transpose:notpossible', ...
    'FOURIETECH objects are not transposable.')
    
end
