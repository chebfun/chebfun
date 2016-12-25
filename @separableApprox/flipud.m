function f = flipud( f ) 
%FLIPUD   Flip/reverse a SEPARABLEAPPROX in the y-direction.
%   G = FLIPUD(F) returns a SEPARABLEAPPROX G with the same domain as F but reversed;
%   that is, G(x,y) = F(x, c+d-y), where the domain is [a, b, c, d].
%
% See also FLIPLR, FLIPDIM. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) )  
    return
end

% Flip the column slices: 
f.cols = flipud( f.cols ); 

end
