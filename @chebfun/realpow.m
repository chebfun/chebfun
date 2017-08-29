function Z = realpow(X, Y)
%REALPOW   Real power of a CHEBFUN.
%   Z = REALPOW(X, Y) is the same as X.^Y.  An error is produced if the result
%   is complex.
%
% See also CHEBFUN/POWER, REALPOW.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

Z = power(X, Y);
        
% Check for complex Z:
if ( ~isreal(Z) )
    error('CHEBFUN:CHEBFUN:realpow:complexRes', ...
        'REALPOW produced complex result.');
end

end
