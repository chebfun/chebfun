function f = fliplr(f) %#ok<INUSD>
%FLIPLR   FLIPLR is not supported for SINGFUN objects.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

error('CHEBFUN:SINGFUN:fliplr:notpossible', ...
      'SINGFUN objects can not be flipped from left to right.')
end
