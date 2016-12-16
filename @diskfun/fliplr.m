function f = fliplr(g)
%FLIPLR   Flip/reverse a DISKFUN over the y-axis. 
%
%   G = FLIPLR( F ) returns a DISKFUN G that is flipped over the y-axis
%   that is, G(x,y) = F(-x,y).
%
% See also DISKFUN/FLIPUD, DISKFUN/FLIPDIM, DISKFUN/ROTATE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( g ) ) 
    return
end 

f = g; 
f.cols = flipud(g.cols); 
f.rows = flipud(g.rows); 

end
