function f = fliplr( f )
%FLIPLR   Flip/reverse a CHEBFUN2 in the x-direction.
%   G = FLIPLR( F ) returns a CHEBFUN2 G with the same domain as F but reversed;
%   that is, G(x,y) = F(a+b-x,y), where the domain is [a, b, c, d].
%
% See also FLIPUD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    return
end

% Flip the row slices: 
f.rows = flipud( f.rows );

end