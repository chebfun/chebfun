function f = imag(f)
%IMAG   Imaginary part of a SINGFUN.
%   IMAG(F) is the imaginary part of F.
%
%   See also REAL, ISREAL, CONJ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Compute the imaginary part of the values:
f.smoothPart = imag(f.smoothPart);

% if F is real, then remove singularities from the imaginary part.
if ( iszero(f.smoothPart) )
    f.exponents = [0, 0];
    f.singType = {'none', 'none'};
end

end
