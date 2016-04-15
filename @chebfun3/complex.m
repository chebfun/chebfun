function h = complex(f, g)
%COMPLEX    Construct complex CHEBFUN3 from real and imaginary parts.
%   H = COMPLEX(F, G) returns the complex CHEBFUN3 F + i G, where F and G 
%   are real valued CHEBFUN3 objects with the same domain.
%
% See also IMAG, CONJ, ABS, REAL.

if ( ~isreal( f ) || ~isreal( g ) )
    error('CHEBFUN:CHEBFUN3:complex:notReal1', ...
        'Inputs must be real valued.');
end
h = f + 1i*g;

end