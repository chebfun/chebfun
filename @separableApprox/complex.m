function C = complex( A, B )
%COMPLEX Construct complex SEPARABLEAPPROX from real and imaginary parts.
%   C = COMPLEX(A, B) returns the complex SEPARABLEAPPROX A + Bi, where A and B are
%   real valued SEPARABLEAPPROX objects with the same domain.
%
%   C = COMPLEX(A) for real SEPARABLEAPPROX A returns the complex result C with real
%   part A and all zero imaginary part. isreal(C) returns false.
%
% See also IMAG, CONJ, ABS, REAL.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 2 )
    if ( ~isa(B, 'separableApprox') )
        error('CHEBFUN:SEPARABLEAPPROX:complex:inputs', ...
            'Second input must be a CHEBFUN2.');
    elseif ( ~isreal( A ) || ~isreal( B ) )
        error('CHEBFUN:SEPARABLEAPPROX:complex:notReal1', ...
            'Inputs must be real valued.');
    end
    C = A + 1i*B;
else
    if ( ~isreal( A ) )
        error('CHEBFUN:SEPARABLEAPPROX:complex:notReal2', ...
            'Input must be real valued.');
    end
    % Make complex.
    C = A + 0*1i;   
end

end
