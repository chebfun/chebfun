function ctranspose(f) %#ok<*INUSD>
%CTRANSPOSE    CLASSICFUN objects are not transposable.
    
% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:CLASSICFUN:ctranspose:notPossible', ...
    'CLASSICFUN objects are not transposable.')
    
end
