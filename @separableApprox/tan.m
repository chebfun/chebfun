function f = tan( f )
%TAN   Tangent of a SEPARABLEAPPROX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    return
end

f = compose( f, @tan ); 

end
