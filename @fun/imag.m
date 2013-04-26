function f = imag(f)
%IMAG	Imaginary part of a FUN.
%   IMAG(F) is the imaginary part of F.
%
%   See also REAL, ISREAL, CONJ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Assign the onefun to be the imaginary part of the onefun of the input:
f.onefun = imag(f.onefun);

end
