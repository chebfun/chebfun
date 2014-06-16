function C = complex( A, B )
%COMPLEX Construct complex CHEBFUN2 from real and imaginary parts.
%   C = COMPLEX(A, B) returns the complex CHEBFUN2 A + Bi, where A and B are
%   real valued CHEBFUN2 objects with the same domain.
%
%   C = COMPLEX(A) for real CHEBFUN2 A returns the complex result C with real
%   part A and all zero imaginary part. isreal(C) returns false.
%
% See also IMAG, CONJ, ABS, REAL.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 2 )
    if ( ~isa(B, 'chebfun2') )
        error('CHEBFUN:CHEBFUN2:complex:inputs', ...
            'Second input must be a CHEBFUN2.');
    elseif ( ~isreal( A ) || ~isreal( B ) )
        error('CHEBFUN:CHEBFUN2:complex:notReal1', ...
            'Inputs must be real valued.');
    end
    C = A + 1i*B;
else
    if ( ~isreal( A ) )
        error('CHEBFUN:CHEBFUN2:complex:notReal2', ...
            'Input must be real valued.');
    end
    % Make complex.
    C = A + 0*1i;   
end

end
