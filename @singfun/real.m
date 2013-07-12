function f = real(f)
%REAL   Real part of a SINGFUN.
%   REAL(F) is the real part of F.
%
%   See also IMAG, ISREAL, CONJ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Compute the imaginary part of the values:
f.smoothPart = real(f.smoothPart);

% if F is imaginary, then remove singularities 
% from the real part
if ( iszero(f.smoothPart) )
    f.exponents = zeros(1,2);
    f.isSingEnd = zeros(1,2);
    f.singType = {'none', 'none'};
end

end
