function f = tanh(f)
%TANH   Hyperbolic tangent of a SEPARABLEAPPROX.
% 
%   TANH(F) returns the hyperbolic tangent of a SEPARABLEAPPROX.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    return
end

f = compose( f, @tanh ); 

end