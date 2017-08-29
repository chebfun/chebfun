function f = sinh( f )
%SINH   Hyperbolic sine of a SEPARABLEAPPROX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) )    
    return
end 

f = compose( f, @sinh ); 

end
