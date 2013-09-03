function f = fliplr(f)
%FLIPLR   Flip columns of an array-valued SINGFUN object.
%   FLIPLR(F) flips the columns of an array-valued SINGFUN in the left/right
%   direction. If F has only one column, then this function has no effect.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

error('CHEBFUN:SINGFUN:fliplr:notpossible', ...
    'SINGFUN objects are not array-valued.')

end
