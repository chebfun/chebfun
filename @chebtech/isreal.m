function out = isreal(f)
%ISREAL   True for real CHEBTECH.
%   ISREAL(F) returns TRUE if F does not have an imaginary part and FALSE
%   otherwise.
%
% See also REAL, IMAG.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if all the coefficients are real:
out = isreal(f.coeffs);

end
