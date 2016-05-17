function f = sin( f ) 
%SIN   Sine of a CHEBFUN2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    return
end

f = compose( f, @sin ); 

end
