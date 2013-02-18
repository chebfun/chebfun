function out = isreal(f)
%ISREAL True for real FUNCHEB2.
%	ISREAL(F) returns 1 if F does not have an imaginary part and 0 otherwise.
%
%   See also REAL, IMAG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Check if all the values are real:
out = isreal(f.values);

end
