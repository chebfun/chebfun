function g = imag(f)
%IMAG   Complex imaginary part of a BALLFUN.
%   IMAG(F) is the imaginary part of F.
%
% See also REAL.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @imag ); 

end
