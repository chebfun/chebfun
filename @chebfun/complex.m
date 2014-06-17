function C = complex(A, B)
%COMPLEX   Construct complex CHEBFUN from real and imaginary parts.
%   COMPLEX(A, B) returns the complex result A + Bi, where A and B are
%   real-valued CHEBFUN objects on the same domain. Alternatively, one of A or
%   B may be a real-valued scalar.
%
% See also REAL, IMAG.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isreal(A) )
    error('CHEBFUN:CHEBFUN:complex:AisNotReal', 'Real input A must be real.');
elseif ( (nargin == 2) && ~isreal(B) )
    error('CHEBFUN:CHEBFUN:complex:BisNotReal', ...
        'Imaginary input B must be real.');
end

if ( nargin == 2 )
    C = A + 1i*B;
else
    C = A;
end

end
