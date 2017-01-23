function f = flipdim(g,dim)
%FLIPDIM   Flip/reverse a DISKFUN in a chosen direction.
%   G = FLIPDIM(F, DIM) returns a DISKFUN G but
%   reversed in a direction:  If DIM = 1 and G(x,y) = F(-x, y).
%   If DIM = 2  then G(x,y) = F(x, -y).  
% 
% See also DISKFUN/FLIPLR, DISKFUN/FLIPUD, DISKFUN/ROTATE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( g ) )  
    return
end

% Use the code in flipud and fliplr. 
if ( dim == 1 ) 
    f = flipud( g ); 
elseif ( dim == 2 )
    f = fliplr( g ); 
else
    error('CHEBFUN:DISKFUN:flipdim:badDim2', 'Dimension not recognised.');
end

end
