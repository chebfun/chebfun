function f = imag(f)
%IMAG   Imaginary part of a CHEBFUN3.
%   IMAG(F) returns the imaginary part of the CHEBFUN3 object F.
%
% See also CHEBFUN3/REAL, CHEBFUN3/CONJ, CHEBFUN3/COMPLEX and 
% CHEBFUN3/COMPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) )
    return
end

f = compose(f, @imag); 

end