function C = rdivide(A, B)
%./    CHEBOP right division.
%   C = A./B, where ones of a A is a CHEBOP and is B a scalar is equivalent
%   to A*(1/B), repectively.
%
%   All other instances returns an error.
%
% See also CHEBOP/MTIMES, CHEBOP/MRDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isnumeric(B) )
    C = mtimes(A, 1./B);
else
    error('CHEBOP:RDIVIDE:NotSupported', '%s./%s is not supported.', ...
        class(A), class(B));
end
    
end
