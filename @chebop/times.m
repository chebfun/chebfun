function C = times(A, B)
%*    CHEBOP multiplication.
%   C = A.*B, where ones of a A or B is a CHEBOP and the other a scalar
%   is equivalent to A*B. 
%
%   All other instances returns an error.
%
% See also CHEBOP/MTIMES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isnumeric(A) || isnumeric(B) )
    C = mtimes(A, B);
else
    error('CHEBOP:TIMES:NotSupported', ...
        ['A.*B is not supported for objects of type A = ' class(A) ...
         ' and B = ' class(B) '.']);
end
    
end
