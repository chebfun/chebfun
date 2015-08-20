function f = imag( f )
%IMAG   Imaginary part of a SEPARABLEAPPROX.
%   IMAG(F) returns the imaginary part of a SEPARABLEAPPROX.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) )
    return
end

f = compose( f, @imag ); 

end