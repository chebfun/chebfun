function f = imag(f)
%IMAG   Imaginary part of a CHEBFUN3.
%   IMAG(F) returns the imaginary part of a CHEBFUN3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) )
    return
end

f = compose(f, @imag); 

end