function C = complex(A, B)
%COMPLEX   Construct complex CHEBFUN from real and imaginary parts.
%   COMPLEX(F, G) returns the complex result F + Gi, where F and G are CHEBFUN
%   objects on the same domain. Alternatively, one of F or G may be a
%   real-valued scalar.
%
% See also REAL, IMAG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isreal(A) )
    error('CHEBFUN:complex:AisNotReal', 'Real input A must be real.');
elseif ( nargin == 2 && ~isreal(B) )
    error('CHEBFUN:complex:BisNotReal', 'Imaginary input B must be real.');
end

if ( nargin == 2 )
    C = A + 1i*B;
else
    C = A;
end

end