function f = tanh(f)
%TANH   Hyperbolic tangent of a CHEBFUN2.
%   TANH(F) returns the hyperbolic tangent of a CHEBFUN2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    return
end

op = @(x,y) tanh( feval( f, x, y ) );  % Resample.    
f = chebfun2( op, f.domain );          % Call constructor.

end