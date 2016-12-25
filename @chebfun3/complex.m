function h = complex(f, g)
%COMPLEX   Construct complex CHEBFUN3 from real and imaginary parts.
%   H = COMPLEX(F, G) returns the complex CHEBFUN3 F + i G, where F and G 
%   are real valued CHEBFUN3 objects with the same domain.
%
% See also CHEBFUN3/IMAG, CHEBFUN3/CONJ, CHEBFUN3/ABS and CHEBFUN3/REAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isreal(f) || ~isreal(g) )
    error('CHEBFUN:CHEBFUN3:complex:notReal1', ...
        'Inputs must be real.');
end
h = f + 1i*g;

end