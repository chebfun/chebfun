function f = sin( f ) 
%SIN   Sine of a CHEBFUN2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    return
end

op = @(x,y) sin( feval(f, x, y) );     % Resample. 
f = chebfun2( op, f.domain );          % Call constructor.

end
