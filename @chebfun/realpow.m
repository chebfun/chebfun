function Z = realpow(X, Y)
%REALSQRT   Real power of a CHEBFUN.
%   Z = REALPOW(X, Y) denotes element-by-element powers. X and Y must have the
%   same dimensions unless one is a scalar. A scalar can operate into anything.
%
%   An error is produced if the result is complex.
%
% See also CHEBFUN/POWER, REALPOW.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

Z = power(X, Y);
        
% Check for complex Z:
if ( ~isreal(Z) )
    error('CHEBFUN:realsqrt:complexR', 'Realpow produced complex result.');
end

end
