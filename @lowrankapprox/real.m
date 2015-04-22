function f = real( f )
%REAL      Real part of a LOWRANKAPPROX.
%
% See also IMAG.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
% Empty check: 
if ( isempty( f ) ) 
    return
end

f = compose( f, @real ); 

end
