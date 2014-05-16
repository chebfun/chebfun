function ctranspose(f) %#ok<*INUSD>
%CTRANSPOSE    CLASSICFUN objects are not transposable.
    
% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:FUN:CLASSICFUN:ctranspose:notpossible', ...
    'CLASSICFUN objects are not transposable.')
    
end