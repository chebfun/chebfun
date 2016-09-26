function f = tand( f )
%TAND  Tangent of a SEPARABLEAPPROX (in degrees)

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) )    
    return
end

f = compose( f, @tand ); 

end
